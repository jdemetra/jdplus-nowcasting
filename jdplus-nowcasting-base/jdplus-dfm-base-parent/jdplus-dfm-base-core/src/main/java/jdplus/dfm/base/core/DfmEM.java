/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jdplus.dfm.base.core;

import java.util.ArrayList;
import jdplus.dfm.base.api.MeasurementType;
import jdplus.dfm.base.api.timeseries.TsInformationSet;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.api.util.Table;
import jdplus.toolkit.base.core.data.DataBlock;

import java.util.EnumMap;
import java.util.List;
import jdplus.dfm.base.core.var.VarDescriptor;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.data.DoubleSeqCursor;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.math.matrices.GeneralMatrix;
import jdplus.toolkit.base.core.math.matrices.SymmetricMatrix;
import jdplus.toolkit.base.core.ssf.ISsfInitialization;
import jdplus.toolkit.base.core.ssf.StateStorage;
import jdplus.toolkit.base.core.ssf.multivariate.MultivariateOrdinaryFilter;
import jdplus.toolkit.base.core.ssf.multivariate.PredictionErrorsDecomposition;
import jdplus.toolkit.base.core.ssf.multivariate.SsfMatrix;
import jdplus.toolkit.base.core.stats.likelihood.Likelihood;

/**
 *
 * @author Jean Palate
 */
public class DfmEM implements IDfmInitializer {

    public static final int MAXITER = 1000, MAXNUMITER = 50;
    public static final double DEF_EPS = 1e-6;

    public static class Builder {

        private IDfmInitializer initializer;
        private int maxIter = MAXITER;
//        private int maxNumericalIter = MAXNUMITER;
        private double eps = DEF_EPS;
//        private boolean fixedInitialConditions = true
        private ISsfInitialization.Type initialization = ISsfInitialization.Type.Zero;
//        private boolean computeAll = true;

        public Builder initializer(IDfmInitializer initializer) {
            this.initializer = initializer;
            return this;
        }

        public Builder maxIter(int maxiter) {
            this.maxIter = maxiter;
            return this;
        }

//        public Builder maxNumericalIter(int maxiter) {
//            this.maxNumericalIter = maxiter;
//            return this;
//        }
//
        public Builder precision(double eps) {
            this.eps = eps;
            return this;
        }

//        public Builder fixedInitialConditions(boolean fic) {
//            this.fixedInitialConditions = fic;
//            return this;
//        }
//
//        public Builder computeAll(boolean all) {
//            this.computeAll = all;
//            return this;
//        }
        public Builder ssfInitialization(ISsfInitialization.Type type) {
            this.initialization = type;
            return this;
        }

        public DfmEM build() {
            return new DfmEM(this);
        }
    }
    
    public static Builder builder(){
        return new Builder();
    }

    private final IDfmInitializer initializer;
    private final int maxIter;
//    private final int maxNumericalIter;
    private final double eps;
//    private final boolean fixedInitialConditions;
//    private final boolean computeAll;
    private final ISsfInitialization.Type initialization;

    private DynamicFactorModel dfm;
    private int nxlags;

    private DfmProcessor processor;
    private TsInformationSet data;
    private Matrix M;
    private final EnumMap<MeasurementType, DataBlock[]> G = new EnumMap<>(MeasurementType.class);
    private final EnumMap<MeasurementType, Table<DataBlock>> G2 = new EnumMap<>(MeasurementType.class);
    private DataBlock Efij[];
    private int iter_;
    private int modelSize;
    private int dataSize;
    private double logLikelihood;

    private DfmEM(Builder builder) {
        this.initializer = builder.initializer;
        this.maxIter = builder.maxIter;
//        this.maxNumericalIter = builder.maxNumericalIter;
        this.eps = builder.eps;
//        this.fixedInitialConditions = builder.fixedInitialConditions;
//        this.computeAll = builder.computeAll;
        this.initialization = builder.initialization;
    }

    public double getFinalLogLikelihood() {
        return logLikelihood;
    }

    /**
     * E(f(i))
     *
     * @param ss
     * @param i
     * @return
     */
    private DataBlock ef(StateStorage ss, int i) {
        return ss.item(i);
    }

    /**
     * cov(f(i), f(j))
     *
     * @param ss
     * @param i
     * @param j
     * @return
     */
    private DataBlock vf(StateStorage ss, int i, int j) {
        return ss.covar(i, j);
    }

    /**
     * Computes E(fi, fj) = cov(fi, fj) + E(fi)*E(fj)
     *
     * @param ss
     * @param i
     * @param j
     * @return
     */
    private DataBlock efifj(StateStorage ss, int i, int j) {
        if (i > j) {
            return efifj(ss, j, i);
        }
        int idx = i + modelSize * j;
        DataBlock cur = Efij[idx];
        if (cur == null) {
            cur = vf(ss, i, j);
            cur.addAXY(1, ef(ss, i), ef(ss, j));
            Efij[idx] = cur;
        }
        return cur;
    }

    /**
     * G, GG contain the smoothed states and their variances aggregated
     * following the different measurement equations
     *
     * @param ss
     */
    private void calcG(StateStorage ss) {
        G.clear();
        G2.clear();
        for (int i = 0; i < Efij.length; ++i) {
            Efij[i] = null;
        }
        int nf = dfm.getNfactors();

        for (MeasurementDescriptor desc : dfm.getMeasurements()) {
            MeasurementType type = IDfmMeasurement.getMeasurementType(desc.getType());
            if (!G.containsKey(type)) {
                int len = desc.getType().getLength();
                DataBlock z = DataBlock.make(len);
                desc.getType().fill(z);
                DataBlock[] g = new DataBlock[nf];
                Table<DataBlock> g2 = new Table<>(nf, nf);
                for (int i = 0, j = 0; i < nf; ++i, j += nxlags) {
                    DataBlock column = DataBlock.make(dataSize);
                    for (int k = 0; k < len; ++k) {
                        column.addAY(z.get(k), ef(ss, j + k));
                    }
                    // g[i] = L
                    g[i] = column;
                }
                for (int i = 0, j = 0; i < nf; ++i, j += nxlags) {
                    for (int k = 0, l = 0; k <= i; ++k, l += nxlags) {
                        DataBlock column = DataBlock.make(dataSize);
                        for (int pr = 0; pr < len; ++pr) {
                            double zr = z.get(pr);
                            if (zr != 0) {
                                for (int pc = 0; pc < len; ++pc) {
                                    double zc = z.get(pc);
                                    if (zc != 0) {
                                        column.addAY(zr * zc, vf(ss, j + pr, l + pc));
                                    }
                                }
                            }
                        }
                        column.addAXY(1, g[i], g[k]);
                        g2.set(i, k, column);
                        if (i != k) {
                            g2.set(k, i, column);
                        }

                    }
                }

                G.put(type, g);
                G2.put(type, g2);
            }
        }
    }

    @Override
    public DynamicFactorModel initialize(DynamicFactorModel rdfm, TsInformationSet data) {
        this.dfm=rdfm;
        this.data = data;
        this.nxlags = rdfm.minSsfLags() + 1;
        this.processor=DfmProcessor.builder()
                .calcVariance(true)
                .initialization(initialization)
                .build();
        modelSize = nxlags * rdfm.getNfactors();
        dataSize = data.getCurrentDomain().getLength();
        Efij = new DataBlock[modelSize * modelSize];
        M = data.generateMatrix(null);
        if (initializer != null) {
            initializer.initialize(dfm, data);
        }
        iter_ = 0;
        logLikelihood = 0;
        filter(true);
        while (iter_++ < maxIter) {
            if (!EStep()) {
                break;
            }
            if (!MStep()) {
                break;
            }
        }

        // finishing
        dfm = dfm.normalize();
        filter(false);
        processor.clear();
        return dfm;
    }

    private void filter(boolean adjust) {
        try {
            MultivariateOrdinaryFilter filter = new MultivariateOrdinaryFilter();
            PredictionErrorsDecomposition results = new PredictionErrorsDecomposition();
            filter.process(dfm.ssfRepresentation(initialization, nxlags), new SsfMatrix(FastMatrix.of(data.generateMatrix(null))), results);
            Likelihood ll = results.likelihood(true);
            logLikelihood = ll.logLikelihood();
            if (adjust) {
                dfm=dfm.rescaleVariances(ll.sigma2());
                dfm=dfm.normalize();
            }
        } catch (RuntimeException err) {
            System.out.println(err.getMessage());
        }
    }

    private boolean EStep() {
        this.processor=DfmProcessor.builder()
                .calcVariance(true)
                .initialization(initialization)
                .extendedLags(nxlags)
                .build();
        if (!processor.process(dfm, data)) {
            return false;
        }
        Likelihood ll = processor.getFilteringResults().likelihood(false);
        if (iter_ > 1 && Math.abs(logLikelihood - ll.logLikelihood()) < eps) {
            return false;
        }
        logLikelihood = ll.logLikelihood();

//        IProcessingHook.HookInformation<DfmEM2, DynamicFactorModel> hinfo
//                = new IProcessingHook.HookInformation<>(this, dfm);
//        this.processHooks(hinfo, all_);
//        if (hinfo.cancel) {
//            return false;
//        }
        calcG(processor.getSmoothingResults());
        return true;
    }

    private boolean MStep() {
        List<MeasurementDescriptor> loading = mloadings();
        VarDescriptor var = mvar();
        DynamicFactorModel tmp = new DynamicFactorModel(var, loading);
        if (tmp.isValid()) {
            dfm = tmp;
            return true;
        } else {
            return false;
        }
    }

    private static double dot(DoubleSeq y, DataBlock g) {
        DoubleSeqCursor ycursor = y.cursor();
        DoubleSeqCursor gcursor = g.cursor();
        int n = y.length();
        double s = 0;
        for (int i = 0; i < n; ++i) {
            double ycur = ycursor.getAndNext();
            if (Double.isFinite(ycur)) {
                s += ycur * gcursor.getAndNext();
            } else {
                gcursor.skip(1);
            }
        }
        return s;
    }

    /**
     * Sum the items of g for non missing y
     *
     * @param y
     * @param g
     * @return
     */
    private static double sum(DoubleSeq y, DataBlock g) {
        DoubleSeqCursor ycursor = y.cursor();
        DoubleSeqCursor gcursor = g.cursor();
        int n = y.length();
        double s = 0;
        for (int i = 0; i < n; ++i) {
            double ycur = ycursor.getAndNext();
            if (Double.isFinite(ycur)) {
                s += gcursor.getAndNext();
            } else {
                gcursor.skip(1);
            }
        }
        return s;
    }

    private List<MeasurementDescriptor> mloadings() {
        // maximize loading
        int i = 0;
        List<MeasurementDescriptor> ndescs = new ArrayList<>();
        // each measurement descriptor corresponds to a variable
        for (MeasurementDescriptor mdesc : dfm.getMeasurements()) {
            MeasurementType type = IDfmMeasurement.getMeasurementType(mdesc.getType());
            DataBlock[] g = G.get(type);
            Table<DataBlock> g2 = G2.get(type);
            // actual variable
            DoubleSeq y = M.column(i++);
            int nobs = y.count(z -> Double.isFinite(z));
            // gy[k] contains E(y*g[k]]
            double[] gy = new double[mdesc.getUsedFactorsCount()];
            FastMatrix g2cur = FastMatrix.square(gy.length);

            for (int j = 0, u = 0; j < mdesc.getCoefficient().length(); ++j) {
                // check that the factor j is used
                if (!Double.isNaN(mdesc.getCoefficient(j))) {
                    // gj
                    DataBlock gj = g[j];
                    gy[u] = dot(y, gj);
                    for (int k = 0, v = 0; k <= j; ++k) {
                        if (!Double.isNaN(mdesc.getCoefficient(k))) {
                            DataBlock g2jk = g2.get(j, k);
                            g2cur.set(u, v, sum(y, g2jk));
                            ++v;
                        }
                    }
                    ++u;
                }
            }
            SymmetricMatrix.fromLower(g2cur);
            // C = gcur*g2cur^-1 or C * g2cur = gy
            SymmetricMatrix.solve(g2cur, DataBlock.of(gy), false);
            double[] c = mdesc.getCoefficient().toArray();
            for (int j = 0, u = 0; j < c.length; ++j) {
                if (!Double.isNaN(c[j])) {
                    c[j] = gy[u++];
                }
            }
            double ee = 0;
            for (int z = 0; z < dataSize; ++z) {
                double yz = y.get(z);
                if (!Double.isNaN(yz)) {
                    ee += yz * yz;
                    for (int j = 0; j < c.length; ++j) {
                        double cj = c[j];
                        if (!Double.isNaN(cj)) {
                            double gj = g[j].get(z);
                            ee -= 2 * cj * gj * yz;
                            for (int k = 0; k < c.length; ++k) {
                                double ck = c[k];
                                if (!Double.isNaN(ck)) {
                                    ee += g2.get(j, k).get(z) * cj * ck;
                                }
                            }
                        }
                    }
                }
            }
            MeasurementDescriptor.Builder mbuilder = mdesc.toBuilder();
            mbuilder.coefficient(DoubleSeq.of(c));
            if (ee < 0) {
                mbuilder.variance(1e-12);
            } else {
                mbuilder.variance(ee / nobs);
            }
            ndescs.add(mbuilder.build());
        }
        return ndescs;
    }

    private VarDescriptor mvar() {
        // analytical optimization
        int nl = dfm.getNlags();
        int nf = dfm.getNfactors();
        int n = nf * nl;
        FastMatrix f = FastMatrix.make(nf, n);
        FastMatrix f2 = FastMatrix.square(n);
        StateStorage ss = processor.getSmoothingResults();
        // fill the matrices
        for (int i = 0; i < nf; ++i) {
            for (int j = 0; j < nl; ++j) {
                for (int k = 0; k < nf; ++k) {
                    double x = efifj(ss, i * nxlags, k * nxlags + j + 1).sum();
                    f.set(i, j * nf + k, x);
                }
            }
        }
        for (int i = 1, r = 0; i <= nl; ++i) {
            for (int k = 0; k < nf; ++k, ++r) {
                for (int j = 1, c = 0; j <= nl; ++j) {
                    for (int l = 0; l < nf; ++l, ++c) {
                        double x = efifj(ss, k * nxlags + i, l * nxlags + j).sum();
                        f2.set(r, c, x);
                    }
                }
            }
        }
        // A = f/f2 <-> Af2 = f
        FastMatrix A = f.deepClone();
        // We don't need f2 anymore ->false
        SymmetricMatrix.solveXS(f2, A, false);

        // Q = 1/T * (E(f0,f0) - A * f')
        FastMatrix Q = FastMatrix.of(dfm.getVar().getInnovationsVariance());
        for (int i = 0; i < nf; ++i) {
            for (int j = 0; j <= i; ++j) {
                Q.set(i, j, efifj(ss, i * nxlags, j * nxlags).sum());
            }
        }
        SymmetricMatrix.fromLower(Q);
        FastMatrix Y = GeneralMatrix.ABt(A, f);
        // Y = A * f'
        Q.sub(Y);
        Q.mul(1.0 / dataSize);

//        if (!fixedInitialConditions && dfm.getInitialization() == ISsfInitialization.Type.Unconditional) {
//            Bfgs bfgs = Bfgs.builder()
//                    .maxIter(maxNumericalIter)
//                    .build();
//            //bfgs.setLineSearch(new SimpleLineSearch());
//            LL2 fn = new LL2();
//            LL2.Instance cur = fn.current();
//            bfgs.minimize(cur);
//            LL2.Instance ofn = (LL2.Instance) bfgs.getResult();
//            if (ofn != null && cur.getValue() > ofn.getValue()) {
//                C = ofn.V;
//                Q = ofn.Q;
//            }
//        }
        return new VarDescriptor(A, Q);
    }

//    // Real function corresponding to the second part of the likelihood
//    class LL2 implements IFunction {
//
//        private final Table<FastMatrix> allK;
//        private final FastMatrix K0;
//
//        LL2() {
//            // computes  results independent of the VAR parameters: K(i,j) = sum(f(i,j), t)
//            VarDescriptor var = dfm.getVarDescriptor();
//            int n = dfm.getBlockLength(), nc = var.getNlags(), nf = dfm.getNfactors();
//            int p = 1 + nc;
//            allK = new Table<>(p, p);
//            for (int i = 0; i < p; ++i) {
//                allK.set(i, i, calcK(i, i));
//                for (int j = 0; j < i; ++j) {
//                    FastMatrix m = calcK(i, j);
//                    allK.set(i, j, m);
//                    allK.set(j, i, GeneralMatrix.transpose(m));
//                }
//            }
////          int n = dfm.getBlockLength(), nc = n-1, nf = dfm.getNfactors();
//            K0 = FastMatrix.square(nc * nf);
//            StateStorage ss = processor.getSmoothingResults();
////           int del = 1;
//            int del = n - nc;
//            for (int i = 0; i < nf; ++i) {
//                for (int k = 0; k < nc; ++k) {
//                    for (int j = 0; j < nf; ++j) {
//                        for (int l = 0; l < nc; ++l) {
//                            double v = efifj(ss, i * n + k + del, j * n + l + del).get(0);
//                            K0.set(i * nc + k, j * nc + l, v);
//                        }
//                    }
//                }
//            }
//        }
//
//        Instance current() {
//            return new Instance();
//        }
//
//        FastMatrix K(int i, int j) {
//            return allK.get(i, j);
//        }
//
//        /**
//         * computes sum(t|f(i+k*len, j+l*len)
//         *
//         * @param i The first lag
//         * @param j The second lag
//         * @return
//         */
//        private FastMatrix calcK(int i, int j) {
//            int n = dfm.getNfactors();
//            int len = dfm.getBlockLength();
//            StateStorage ss = processor.getSmoothingResults();
//            FastMatrix K = FastMatrix.square(n);
//            int nlags = dfm.getNlags();
//            for (int k = 0; k < n; ++k) {
//                for (int l = 0; l < n; ++l) {
//                    double s = efifj(ss, i + k * len, j + l * len).sum();
////                    // add first ef...
//                    for (int u = 1; u < len - nlags; ++u) {
//                        s += efifj(ss, i + u + k * len, j + u + l * len).get(0);
//                    }
//                    K.set(k, l, s);
//                }
//            }
//            return K;
//        }
//
//        @Override
//        public Instance evaluate(DoubleSeq parameters) {
//            return new Instance(parameters);
//        }
//
//        @Override
//        public IParametersDomain getDomain() {
//            int nf = dfm.getNfactors();
//            int nl = dfm.getNlags();
//            return new DefaultDomain(nf * nf * nl + nf * (nf + 1) / 2, 1e-6);
//        }
//
//        public class Instance implements IFunctionPoint {
//
//            private final FastMatrix Q, lQ, V; // lQ = cholesky factor of Q
//            private final DataBlock p;
//            private final double val;
//            private final FastMatrix lv0;
//            private final FastMatrix[] LA;
//            private double v0, v1, v2, v3;
//
//            public Instance(DoubleSeq p) {
//                this.p = DataBlock.of(p);
//                VarDescriptor var = dfm.getVar();
//                int nf = var.getNfactors(), nl = var.getNlags();
//                V = FastMatrix.make(nf, nf * nl);
//                Q = FastMatrix.make(nf, nf * nl);
//                int vlen = V.size();
//                this.p.range(0, vlen).copyTo(V.getStorage(), 0);
//                for (int i = 0, j = vlen; i < nf; j += nf - i, i++) {
//                    Q.column(i).drop(i, 0).copy(this.p.range(j, j + nf - i));
//                }
//                SymmetricMatrix.fromLower(Q);
//                lQ = Q.deepClone();
//                SymmetricMatrix.lcholesky(lQ);
//                lv0 = calclv0();
//                LA = calcLA();
//                val = calc();
//            }
//
//            public Instance() {
//                @NonNull VarDescriptor var = dfm.getVar();
//                this.V = FastMatrix.of(var.getCoefficients());
//                this.Q = FastMatrix.of(var.getInnovationsVariance());
//                int nf = dfm.getNfactors();
//                int vlen = V.size();
//                double[] buffer = new double[vlen + nf * (nf + 1) / 2];
//                V.copyTo(buffer, 0);
//                p = DataBlock.of(buffer);
//                for (int i = 0, j = vlen; i < nf; j += nf - i, i++) {
//                    p.range(j, j + nf - i).copy(Q.column(i).drop(i, 0));
//                }
//                this.lQ = Q.deepClone();
//                SymmetricMatrix.lcholesky(lQ);
//                lv0 = calclv0();
//                LA = calcLA();
//                val = calc();
//            }
//
//            private double calc() {
//                v0 = calcdetv0();
//                v1 = calcssq0();
//                v2 = calcdetq();
//                v3 = calcssq();
//                return v0 + v1 + v2 + v3;
//            }
//
//            private FastMatrix A(int i) {
//                VarDescriptor var = dfm.getVar();
//                int nf = var.getNfactors(), nl = var.getNlags();
//                if (i == 0) {
//                    return FastMatrix.identity(nf);
//                } else {
//                    FastMatrix a = FastMatrix.square(nf);
//                    for (int j = 0; j < nf; ++j) {
//                        a.column(j).copy(V.column(i - 1 + j * nl));
//                    }
//                    a.chs();
//                    return a;
//                }
//            }
//
//            @Override
//            public IFunctionDerivatives derivatives() {
//                return new NumericalDerivatives(this, false, true);
//            }
//
//            @Override
//            public DoubleSeq getParameters() {
//                return p;
//            }
//
//            @Override
//            public double getValue() {
//                return val; //To change body of generated methods, choose Tools | Templates.
//            }
//
//            private double calcssq0() {
//                // computes f0*f0 x V^-1
//                // V^-1 = (LL')^-1 = L'^-1*L^-1
//                FastMatrix lower = LowerTriangularMatrix.inverse(lv0);
//                FastMatrix iv0 = SymmetricMatrix.LtL(lower);
//                return iv0.dot(K0);
//
//            }
//
//            private double calcdetv0() {
//                LogSign sumLog = LogSign.of(lv0.diagonal());
//                if (!sumLog.isPositive()) {
//                    throw new DfmException();
//                }
//                return 2 * sumLog.getValue();
//            }
//
//            private double calcdetq() {
//                LogSign sumLog = LogSign.of(lQ.diagonal());
//                if (!sumLog.isPositive()) {
//                    throw new DfmException();
//                }
//                return (dataSize + dfm.getBlockLength() - dfm.getNlags() - 1) * 2 * sumLog.getValue();
////               return dataSize * 2 * sumLog.value;
//            }
//
//            private double calcssq() {
//                int nl = 1 + dfm.getVarDescriptor().getNlags();
//                int nf = dfm.getNfactors();
//                double ssq = 0;
//                FastMatrix aqa = FastMatrix.square(nf);
//                for (int i = 0; i < nl; ++i) {
//                    SymmetricMatrix.XtX(LA[i], aqa);
//                    ssq += K(i, i).dot(aqa);
//                    for (int j = 0; j < i; ++j) {
//                        GeneralMatrix.setAtB(LA[i], LA[j], aqa);
//                        ssq += 2 * K(i, j).dot(aqa);
//                    }
//                }
//                return ssq;
//            }
//
//            private FastMatrix calclv0() {
//                switch (dfm.getInitialization()) {
//                    case UserDefined -> {
//                        return FastMatrix.of(dfm.getInitialVariance());
//                    }
//                    case Unconditional -> {
//                        // compute the initial covar. We reuse the code of DynamicFactorModel
//                        DynamicFactorModel.Builder builder = dfm.toBuilder();
//                        DynamicFactorModel tmp = builder.clearMeasurements()
//                                .nlags(dfm.getVarDescriptor().getNlags())
//                                .var(dfm.getVar().toBuilder()
//                                        .coefficients(V)
//                                        .innovationsVariance(Q)
//                                        .buildWithoutValidation())
//                                .build();
//
//                        try {
//                            int n = tmp.getNfactors() * tmp.getBlockLength();
//                            FastMatrix cov = FastMatrix.square(n);
//                            tmp.ssfRepresentation().initialization().Pf0(cov);
//                            SymmetricMatrix.lcholesky(cov);
//                            LowerTriangularMatrix.toLower(cov);
//                            return cov;
//                        } catch (MatrixException err) {
//                            throw new DfmException();
//                        }
//                    }
//                    default -> {
//                        return null;
//                    }
//                }
//
//            }
//
//            /**
//             * AQA = A'(i)*Q^-1*A(j) = A"(i)*L'^-1*L^-1*A We compute here L^-1*A
//             */
//            private FastMatrix[] calcLA() {
//                FastMatrix[] M = new FastMatrix[1 + dfm.getVar().getNlags()];
//                for (int i = 0; i < M.length; ++i) {
//                    FastMatrix ai = A(i);
//                    // L^-1*A = B or  L B = A
//                    LowerTriangularMatrix.solveLX(lQ, ai);
//                    M[i] = ai;
//                }
//                return M;
//            }
//
//            @Override
//            public IFunction getFunction() {
//                return LL2.this;
//            }
//        }
//    }
//
}

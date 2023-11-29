/*
 * Copyright 2023 National Bank of Belgium
 * 
 * Licensed under the EUPL, Version 1.2 or – as soon they will be approved 
 * by the European Commission - subsequent versions of the EUPL (the "Licence");
 * You may not use this work except in compliance with the Licence.
 * You may obtain a copy of the Licence at:
 * 
 * https://joinup.ec.europa.eu/software/page/eupl
 * 
 * Unless required by applicable law or agreed to in writing, software 
 * distributed under the Licence is distributed on an "AS IS" basis,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the Licence for the specific language governing permissions and 
 * limitations under the Licence.
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
import jdplus.dfm.base.api.EmSpec;
import jdplus.dfm.base.core.var.VarDescriptor;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.data.DoubleSeqCursor;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.math.matrices.GeneralMatrix;
import jdplus.toolkit.base.core.math.matrices.SymmetricMatrix;
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
//        private boolean computeAll = true;
        private boolean fixedVar = false;

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

        public Builder fixedVar(boolean fixed) {
            this.fixedVar = fixed;
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
        public DfmEM build() {
            return new DfmEM(this);
        }
    }

    public static Builder builder() {
        return new Builder();
    }

    public static DfmEM of(EmSpec spec) {
        return builder()
                .maxIter(spec.getMaxIter())
                .precision(spec.getPrecision())
                .build();
    }

    private final IDfmInitializer initializer;
    private final int maxIter;
//    private final int maxNumericalIter;
    private final double eps;
//    private final boolean fixedInitialConditions;
//    private final boolean computeAll;
    private final boolean fixedVar;

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
        this.fixedVar = builder.fixedVar;
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
        this.dfm = rdfm;
        this.data = data;
        this.nxlags = rdfm.minSsfBlockLength();
        this.processor = DfmProcessor.builder()
                .calcVariance(true)
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
            filter.process(dfm.ssfRepresentation(0), new SsfMatrix(FastMatrix.of(data.generateMatrix(null))), results);
            Likelihood ll = results.likelihood(true);
            logLikelihood = ll.logLikelihood();
            if (adjust) {
                dfm = dfm.rescaleVariances(ll.sigma2());
                dfm = dfm.normalize();
            }
        } catch (RuntimeException err) {
        }
    }

    private boolean EStep() {
        this.processor = DfmProcessor.builder()
                .calcVariance(true)
                .extendedLags(0)
                .build();
        if (!processor.process(dfm, data)) {
            return false;
        }
        Likelihood ll = processor.getFilteringResults().likelihood(true);
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
        if (fixedVar || tmp.isValid()) {
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
        VarDescriptor var = dfm.getVar();
        if (fixedVar) {
            return var;
        }
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
        FastMatrix Q = FastMatrix.of(var.getInnovationsVariance()
        );
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

        return new VarDescriptor(A, Q, var.getInitialization());
    }

}

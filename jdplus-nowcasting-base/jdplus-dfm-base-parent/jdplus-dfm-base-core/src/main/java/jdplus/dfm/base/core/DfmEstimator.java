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

import internal.jdplus.dfm.base.core.DfmFunction;
import internal.jdplus.dfm.base.core.DfmFunctionPoint;
import jdplus.dfm.base.api.NumericalProcessingSpec;
import jdplus.dfm.base.api.timeseries.TsInformationSet;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.data.DoublesMath;
import jdplus.toolkit.base.api.information.GenericExplorable;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.api.timeseries.TsDomain;
import jdplus.toolkit.base.core.data.DataBlock;
import jdplus.toolkit.base.core.data.DataBlockIterator;
import jdplus.toolkit.base.core.math.functions.FunctionMinimizer;
import jdplus.toolkit.base.core.math.functions.bfgs.Bfgs;
import jdplus.toolkit.base.core.math.functions.levmar.LevenbergMarquardtMinimizer;
import jdplus.toolkit.base.core.math.functions.ssq.ProxyMinimizer;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.ssf.ISsfInitialization;
import jdplus.toolkit.base.core.ssf.multivariate.SsfMatrix;
import jdplus.toolkit.base.core.stats.likelihood.Likelihood;

/**
 *
 * @author Jean Palate
 */
public class DfmEstimator implements IDfmEstimator {

    private static final String SIMPLIFIED = "Optimizing simplified model",
            MSTEP = "Optimizing measurements", VSTEP = "Optimizing Var model", ALL = "Optimizing all parameters";

    public static final int DEF_MAXITER = 100, DEF_NSTART = 0, DEF_NNEXT = 5, DEF_MAXBLOCKITERATIONS = 100, DEF_MAXEMUP = 5;

    public static class Builder {

        private int maxIter = DEF_MAXITER, maxBlockIterations = DEF_MAXBLOCKITERATIONS, maxEmUp = DEF_MAXEMUP;
        private FunctionMinimizer.Builder minimizer;
        private int maxInitialIter = DEF_NSTART, maxIntermediateIter = DEF_NNEXT;
        private boolean mixed = true, independentVarShocks = true;
        private double eps = 1e-9;
        private TsDomain edomain = null;

        public Builder maxIterations(int maxiter) {
            this.maxIter = maxiter;
            return this;
        }

        public Builder maxInitialIter(int maxInitialIter) {
            this.maxInitialIter = maxInitialIter;
            return this;
        }

        public Builder maxIntermediateIter(int maxIntermediateIter) {
            this.maxIntermediateIter = maxIntermediateIter;
            return this;
        }

        public Builder maxEmUp(int maxEmUp) {
            this.maxEmUp = maxEmUp;
            return this;
        }

        public Builder maxBlockIterations(int maxblocks) {
            this.maxBlockIterations = maxblocks;
            return this;
        }

        public Builder mixed(boolean mixed) {
            this.mixed = mixed;
            return this;
        }

        public Builder independentVarShocks(boolean independentVarShocks) {
            this.independentVarShocks = independentVarShocks;
            return this;
        }

        public Builder minimizer(FunctionMinimizer.Builder minimizer) {
            this.minimizer = minimizer;
            return this;
        }

        public Builder estimationDomain(TsDomain domain) {
            this.edomain = domain;
            return this;
        }

        public Builder precision(double eps) {
            this.eps = eps;
            return this;
        }

        public DfmEstimator build() {
            return new DfmEstimator(this);
        }
    }

    public static Builder builder() {
        return new Builder();
    }

    public static DfmEstimator of(NumericalProcessingSpec spec, ISsfInitialization.Type type) {

        FunctionMinimizer.Builder minimizer;
        minimizer = switch (spec.getMethod()) {
            case BFGS ->
                Bfgs.builder();
            default ->
                ProxyMinimizer.builder(LevenbergMarquardtMinimizer.builder());
        };

        return builder()
                .minimizer(minimizer)
                .maxInitialIter(spec.getMaxInitialIter())
                .maxIntermediateIter(spec.getMaxIntermediateIter())
                .mixed(spec.isMixedEstimation())
                .maxBlockIterations(spec.isEstimationByBlock() ? DEF_MAXBLOCKITERATIONS : 0)
                .independentVarShocks(spec.isIndependentShocks())
                .precision(spec.getPrecision())
                .build();
    }

    private final int maxIter;
    private final int maxInitialIter, maxIntermediateIter, maxBlockIterations, maxEmUp;
    private final boolean mixed;
    private final boolean independentVarShocks;
    private final double eps;
    private final FunctionMinimizer.Builder minimizer;
    private final TsDomain edomain;

    private boolean converged;
    private DynamicFactorModel dfm;
    private Likelihood likelihood;
    private FastMatrix hessian;
    private DoubleSeq gradient;

    public DfmEstimator(Builder builder) {
        this.maxIter = builder.maxIter;
        this.maxInitialIter = builder.maxInitialIter;
        this.maxIntermediateIter = builder.maxIntermediateIter;
        this.maxBlockIterations = builder.maxBlockIterations;
        this.maxEmUp = builder.maxEmUp;
        this.independentVarShocks = builder.independentVarShocks;
        this.mixed = builder.mixed;
        this.minimizer = builder.minimizer;
        this.edomain = builder.edomain;
        this.eps = builder.eps;
    }

    public Builder toBuilder() {
        return builder()
                .maxIterations(maxIter)
                .maxInitialIter(maxInitialIter)
                .maxIntermediateIter(maxIntermediateIter)
                .maxBlockIterations(maxBlockIterations)
                .maxEmUp(maxEmUp)
                .estimationDomain(edomain)
                .independentVarShocks(independentVarShocks)
                .mixed(mixed)
                .minimizer(minimizer)
                .precision(eps);

    }

    public boolean hasConverged() {
        return converged;
    }

    private DynamicFactorModel normalize(DynamicFactorModel model) {
        if (independentVarShocks) {
            return model.lnormalize();
        } else {
            return model.normalize();
        }
    }

    private IDfmMapping mapping(DynamicFactorModel model, boolean mf, boolean vf) {
        if (independentVarShocks) {
            return new DfmMappingI(model, mf, vf);
        } else {
            return new DfmMapping(model, mf, vf);
        }
    }

    @Override
    public boolean estimate(final DynamicFactorModel dfm, TsInformationSet input) {
        DynamicFactorModel model = dfm;
        converged = false;
        SsfMatrix data = new SsfMatrix(FastMatrix.of(input.generateMatrix(edomain)));
        DfmFunction fn;
        DfmFunctionPoint pt;
        model = model.normalize();
        FunctionMinimizer fnmin = null;
        boolean log = false;
        int emUpLeft = maxEmUp;
        int niter = 0;
        try {
            if (maxInitialIter > 0) {
//                setMessage(SIMPLIFIED);
                fnmin = minimizer
                        .maxIter(maxInitialIter)
                        .functionPrecision(eps)
                        .build();
                SimpleDfmMapping smapping = new SimpleDfmMapping(model);
                DynamicFactorModel smodel = smapping.validate(model);

                fn = DfmFunction.builder(data, smapping)
                        .parallelProcessing(true)
                        .symmetricNumericalDerivatives(false)
                        .log(log)
                        .build();
                DfmFunctionPoint curpt = fn.evaluate(smapping.map(smodel));
//                System.out.println(curpt.getLikelihood().logLikelihood());
                fnmin.minimize(curpt);
                pt = (DfmFunctionPoint) fnmin.getResult();
                likelihood = pt.getLikelihood();
//                System.out.println(pt.getLikelihood().logLikelihood());
                double var = likelihood.sigma2();
                model = pt.getCore().rescaleVariances(var);
            }
            if (maxBlockIterations > 0) {
                fnmin = minimizer
                        .maxIter(maxIntermediateIter)
                        .functionPrecision(eps)
                        .build();
                while (true) {
                    model = normalize(model);
                    IDfmMapping mapping = mapping(model, true, false);
                    fn = DfmFunction.builder(data, mapping)
                            .parallelProcessing(true)
                            .symmetricNumericalDerivatives(false)
                            .log(log)
                            .build();

//                    setMessage(VSTEP);
                    pt = fn.evaluate(mapping.map(model));
                    fnmin.minimize(pt);
                    niter += fnmin.getIterationsCount();
                    pt = (DfmFunctionPoint) fnmin.getResult();
                    likelihood = pt.getLikelihood();

//                    System.out.println(pt.getLikelihood().logLikelihood());
                    double var = likelihood.sigma2();
                    model = pt.getCore().rescaleVariances(var);
                    model = normalize(model);
                    if (mixed) {
                        double ll0 = pt.getLikelihood().logLikelihood();
                        DfmEM em = DfmEM.builder()
                                .maxIter(maxIntermediateIter)
                                .fixedVar(emUpLeft <= 0)
                                .build();
                        model = em.initialize(model, input);
                        double ll1 = em.getFinalLogLikelihood();
                        if (ll1 < ll0) {
                            --emUpLeft;
                        }
                    } else {
                        mapping = mapping(model, false, true);
                        fn = DfmFunction.builder(data, mapping)
                                .parallelProcessing(true)
                                .symmetricNumericalDerivatives(false)
                                .log(log)
                                .build();
//                        setMessage(MSTEP);
                        fnmin.minimize(fn.evaluate(mapping.map(model)));
                        niter += fnmin.getIterationsCount();
                        pt = (DfmFunctionPoint) fnmin.getResult();
//                        System.out.println(pt.getLikelihood().logLikelihood());
                        var = pt.getLikelihood().sigma2();
                        model = pt.getCore().rescaleVariances(var);
                        model = normalize(model);
                    }
                    mapping = mapping(model, false, false);
                    fn = DfmFunction.builder(data, mapping)
                            .parallelProcessing(true)
                            .symmetricNumericalDerivatives(false)
                            .log(log)
                            .build();
//                    setMessage(ALL);
                    converged = fnmin.minimize(fn.evaluate(mapping.map(model)));
                    niter += fnmin.getIterationsCount();
                    pt = (DfmFunctionPoint) fnmin.getResult();
//                    System.out.println(pt.getLikelihood().logLikelihood());
                    var = pt.getLikelihood().sigma2();
                    model = pt.getCore().rescaleVariances(var);
                    model = normalize(model);
                    boolean stop = likelihood != null && Math.abs(likelihood.logLikelihood() - pt.getLikelihood().logLikelihood()) < eps;
                    likelihood = pt.getLikelihood();
                    if (converged || niter >= maxIter || stop) {
                        break;
                    }
                }
            } else {
                model = normalize(model);
                IDfmMapping mapping = mapping(model, false, false);
                fn = DfmFunction.builder(data, mapping)
                        .parallelProcessing(true)
                        .symmetricNumericalDerivatives(false)
                        .log(log)
                        .build();
                fnmin = minimizer
                        .maxIter(maxIter)
                        .functionPrecision(eps)
                        .build();
//                setMessage(ALL);
                converged = fnmin.minimize(fn.evaluate(mapping.map(model)));
                pt = (DfmFunctionPoint) fnmin.getResult();
                double var = pt.getLikelihood().sigma2();
                model = pt.getCore().rescaleVariances(var);
                likelihood = pt.getLikelihood();
            }
            return true;
        } catch (Exception err) {
            return false;
        } finally {
            model = normalize(model);
            this.dfm = model;
            if (fnmin != null) {
                finalizeProcessing(fnmin);
            }
        }
    }

    private void finalizeProcessing(FunctionMinimizer fnmin) {
        IDfmMapping fmapping = mapping(dfm, false, false);
        DoubleSeq mp = fmapping.map(dfm);
        DoubleSeq up = fnmin.getResult().getParameters();
        DoubleSeq factors = DoublesMath.divide(mp, up);

        DfmFunction function = (DfmFunction) fnmin.getResult().getFunction();
        boolean log = function.isLog();

        FastMatrix h = fnmin.curvatureAtMinimum();
        if (h != null) {
            if (!log) {
                // we have to correct the hessian 
                int ndf = likelihood.dim() - h.getRowsCount();
                h = h.times(.5 * ndf / fnmin.getObjective());
                DataBlockIterator rows = h.columnsIterator(), cols = h.columnsIterator();
                while (rows.hasNext()) {
                    rows.next().mul(factors);
                }
                while (cols.hasNext()) {
                    cols.next().mul(factors);
                }
            }
            hessian = h;
        }
        DataBlock grad = DataBlock.of(fnmin.gradientAtMinimum());
        if (grad != null) {
            if (!log) {
                // we have to correct the gradient 
                int ndf = likelihood.dim() - grad.length();
                grad.mul(-.5 * ndf / fnmin.getObjective());
                grad.mul(factors);
            }
            gradient = grad;
        }
    }

    @Override
    public Matrix getHessian() {
        return this.hessian;
    }

    @Override
    public DoubleSeq getGradient() {
        return this.gradient;
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public DynamicFactorModel getEstimatedModel() {
        return dfm;
    }
}

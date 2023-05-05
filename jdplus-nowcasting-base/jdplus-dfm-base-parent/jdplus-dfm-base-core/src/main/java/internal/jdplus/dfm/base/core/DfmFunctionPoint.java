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
package internal.jdplus.dfm.base.core;

import jdplus.dfm.base.core.DynamicFactorModel;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.core.math.functions.IFunction;
import jdplus.toolkit.base.core.math.functions.IFunctionDerivatives;
import jdplus.toolkit.base.core.math.functions.NumericalDerivatives;
import jdplus.toolkit.base.core.math.functions.ssq.ISsqFunction;
import jdplus.toolkit.base.core.math.functions.ssq.ISsqFunctionDerivatives;
import jdplus.toolkit.base.core.math.functions.ssq.SsqNumericalDerivatives;
import jdplus.toolkit.base.core.ssf.SsfException;
import jdplus.toolkit.base.core.ssf.multivariate.IMultivariateSsf;
import jdplus.toolkit.base.core.ssf.multivariate.MultivariateOrdinaryFilter;
import jdplus.toolkit.base.core.stats.likelihood.Likelihood;
import jdplus.toolkit.base.core.stats.likelihood.LikelihoodFunctionPoint;

/**
 *
 * @author palatej
 */
public class DfmFunctionPoint implements
        LikelihoodFunctionPoint<Likelihood> {

    /**
     *
     */
    private final IMultivariateSsf currentSsf;
    private final DynamicFactorModel current;

    /**
     *
     */
    private final Likelihood ll;
    private final DoubleSeq p;
    private final DfmFunction fn;

    /**
     *
     * @param fn
     * @param p
     */
    public DfmFunctionPoint(DfmFunction fn, DoubleSeq p) {
        this.fn = fn;
        this.p = p;
        current = fn.getMapping().map(p);
        currentSsf = current.ssfRepresentation(0);
        Likelihood l=null;
        try {
            MultivariateOrdinaryFilter filter = new MultivariateOrdinaryFilter();
            PredictionErrorsDecompositionEx results = new PredictionErrorsDecompositionEx();
            filter.process(currentSsf, fn.getData(), results);
            l = results.likelihood(true);
        } catch (SsfException err) {
        }
        ll=l;
    }

    public DynamicFactorModel getCore() {
        return current;
    }

    @Override
    public DoubleSeq getE() {
        return ll == null ? null : ll.deviances();
    }

    /**
     *
     * @return
     */
    @Override
    public Likelihood getLikelihood() {
        return ll;
    }

    @Override
    public DoubleSeq getParameters() {
        return p;
    }

    @Override
    public double getSsqE() {
        if (ll == null) {
            return Double.NaN;
        }
        return ll.ssq() * ll.factor();
    }

    @Override
    public double getValue() {
        if (ll == null) {
            return Double.NaN;
        }
        if (fn.isLog()) {
            return -ll.logLikelihood();
        } else {
            return ll.ssq() * ll.factor();
        }
    }

    @Override
    public ISsqFunction getSsqFunction() {
        return fn;
    }

    @Override
    public IFunction getFunction() {
        return fn;
    }

    @Override
    public IFunctionDerivatives derivatives() {
        return new NumericalDerivatives(this, fn.isSymmetric(), fn.isMultiThreaded());
    }

    @Override
    public ISsqFunctionDerivatives ssqDerivatives() {
        return new SsqNumericalDerivatives(this, fn.isSymmetric(), fn.isMultiThreaded());
    }
}

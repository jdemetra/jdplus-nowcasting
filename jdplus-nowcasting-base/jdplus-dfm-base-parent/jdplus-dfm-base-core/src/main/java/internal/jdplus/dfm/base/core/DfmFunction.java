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

import jdplus.dfm.base.core.IDfmMapping;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.core.math.functions.IParametersDomain;
import jdplus.toolkit.base.core.ssf.ISsfInitialization;
import jdplus.toolkit.base.core.ssf.multivariate.IMultivariateSsfData;
import jdplus.toolkit.base.core.stats.likelihood.Likelihood;
import jdplus.toolkit.base.core.stats.likelihood.LikelihoodFunction;
import nbbrd.design.BuilderPattern;

/**
 *
 * @author palatej
 */
public class DfmFunction implements LikelihoodFunction<Likelihood> {

    @BuilderPattern(DfmFunction.class)
    public static class Builder {

        private final IDfmMapping mapping;
        private final IMultivariateSsfData data;
        private boolean log = false, mt = false, sym = false;
        private ISsfInitialization.Type initialization = ISsfInitialization.Type.Unconditional;

        private Builder(final IMultivariateSsfData data, final IDfmMapping mapping) {
            this.data = data;
            this.mapping = mapping;
        }

        public Builder parallelProcessing(boolean mt) {
            this.mt = mt;
            return this;
        }

        public Builder log(boolean log) {
            this.log = log;
            return this;
        }

        public Builder initialization(ISsfInitialization.Type initialization) {
            this.initialization = initialization;
            return this;
        }

        public Builder symmetricNumericalDerivatives(boolean sym) {
            this.sym = sym;
            return this;
        }

        public DfmFunction build() {
            return new DfmFunction(this);
        }
    }

    public static Builder builder(IMultivariateSsfData data, IDfmMapping mapping) {
        return new Builder(data, mapping);
    }

    private final IDfmMapping mapping; // mapping from an array of double to an object S
    private final IMultivariateSsfData data;
    private final boolean log, mt, sym;
    private final ISsfInitialization.Type initialization;

    private DfmFunction(Builder builder) {
        this.data = builder.data;
        this.mapping = builder.mapping;
        this.log = builder.log;
        this.mt = builder.mt;
        this.sym = builder.sym;
        this.initialization = builder.initialization;
    }

    @Override
    public DfmFunctionPoint evaluate(DoubleSeq parameters) {
        return new DfmFunctionPoint(this, parameters);
    }

    /**
     *
     * @return
     */
    @Override
    public IParametersDomain getDomain() {
        return mapping;
    }

    @Override
    public DfmFunctionPoint ssqEvaluate(DoubleSeq parameters) {
        return new DfmFunctionPoint(this, parameters);
    }

    /**
     * @return the data
     */
    public IMultivariateSsfData getData() {
        return data;
    }

    public IDfmMapping getMapping() {
        return mapping;
    }

    public ISsfInitialization.Type getInitialization() {
        return initialization;
    }

    /**
     * @return the mt
     */
    public boolean isMultiThreaded() {
        return mt;
    }

    /**
     * @return the sym
     */
    public boolean isSymmetric() {
        return sym;
    }

    public boolean isLog() {
        return log;
    }

}

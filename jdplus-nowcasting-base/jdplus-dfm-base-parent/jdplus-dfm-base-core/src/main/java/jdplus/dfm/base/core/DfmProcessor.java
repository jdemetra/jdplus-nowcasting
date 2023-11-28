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

import internal.jdplus.dfm.base.core.MSmoother;
import jdplus.dfm.base.api.DfmException;
import jdplus.dfm.base.api.timeseries.TsInformationSet;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.ssf.StateStorage;
import jdplus.toolkit.base.core.ssf.multivariate.IMultivariateSsf;
import jdplus.toolkit.base.core.ssf.multivariate.MultivariateFilteringInformation;
import jdplus.toolkit.base.core.ssf.multivariate.SsfMatrix;

/**
 *
 * @author Jean Palate
 */
public class DfmProcessor {

    public static class Builder {

        private boolean calcVariance = true;
        private int nxlags;

        public Builder calcVariance(boolean var) {
            this.calcVariance = var;
            return this;
        }

        public Builder extendedLags(int nxlags) {
            this.nxlags = nxlags;
            return this;
        }

        public DfmProcessor build(){
            return new DfmProcessor(this);
        }
        
    }
    
    public static Builder builder(){
        return new Builder();
    }

    private final boolean calcVariance;
    private final int nxlags;
    
    private DfmProcessor(Builder builder){
        this.calcVariance=builder.calcVariance;
        this.nxlags=builder.nxlags;
    }

    private StateStorage smoothingResults;
    private MultivariateFilteringInformation filteringResults;
    public void clear() {
        smoothingResults = null;
        filteringResults = null;
    }

     /**
     * Retrieves the smoothing results
     *
     * @return The Smoothing results. May by null.
     */
    public StateStorage getSmoothingResults() {
        return smoothingResults;
    }

    public MultivariateFilteringInformation getFilteringResults() {
        return filteringResults;
    }

    public boolean process(DynamicFactorModel model, TsInformationSet input) {
        Matrix M = input.generateMatrix(null);
        if (M.getColumnsCount() != model.getMeasurementsCount()) {
            throw new DfmException(DfmException.INCOMPATIBLE_DATA);
        }
        try {
            clear();
            IMultivariateSsf ssf = model.ssfRepresentation(nxlags);
            MSmoother smoother = MSmoother.builder(ssf)
                    .calcVariance(calcVariance)
                    .build();

            if (smoother.process(new SsfMatrix(FastMatrix.of(M)))) {
                smoothingResults = smoother.getSmoothingResults();
                filteringResults = smoother.getFilteringResults();
                return true;
            } else {
                return false;
            }
        } catch (RuntimeException err) {
            smoothingResults = null;
            return false;
        }
    }

}

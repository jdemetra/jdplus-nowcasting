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

import jdplus.dfm.base.api.timeseries.TsInformationSet;
import jdplus.toolkit.base.core.ssf.StateStorage;
import jdplus.toolkit.base.core.ssf.multivariate.MultivariateFilteringInformation;

/**
 *
 * @author Jean Palate
 */
@lombok.Value
@lombok.Builder(builderClassName="Builder")
public class DfmKernel {

    private IDfmInitializer initializer;
    private IDfmEstimator estimator;
    private DfmProcessor processor;

    public DfmResults process(DynamicFactorModel model, TsInformationSet input) {
        DfmResults.builder builder = DfmResults.builder();
        DynamicFactorModel dfm=model;
        if (initializer != null) {
            dfm=initializer.initialize(dfm, input);
            builder.dfm(dfm);
        }
        if (estimator != null) {
            if (!estimator.estimate(dfm, input)) {
                return null;
            }
            dfm=estimator.getEstimatedModel();
            builder.dfm(dfm)
                    .hessian(estimator.getHessian())
                    .gradient(estimator.getGradient());
        }
        if (processor != null) {
            processor.process(dfm, input);
            builder.smoothedStates(processor.getSmoothingResults());
        }
        
        return builder.build();
    }

    public StateStorage getSmoothingResults() {
        return processor != null ? processor.getSmoothingResults() : null;
    }
    public MultivariateFilteringInformation getFilteringResults() {
        return processor != null ? processor.getFilteringResults(): null;
    }
}

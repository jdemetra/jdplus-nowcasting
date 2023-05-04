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

import jdplus.dfm.base.api.NumericalProcessingSpec;
import static jdplus.dfm.base.core.DfmEMTest.dfmdata;
import static jdplus.dfm.base.core.DfmEMTest.dmodel;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.timeseries.TsDomain;
import jdplus.toolkit.base.core.ssf.ISsfInitialization;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

/**
 *
 * @author palatej
 */
public class DfmEstimatorTest {

    public DfmEstimatorTest() {
    }

    public static void main(String[] args) {

        DfmEM em = DfmEM.builder()
                .ssfInitialization(ISsfInitialization.Type.Unconditional)
                .maxIter(10)
                .build();
        PrincipalComponentsInitializer initializer = new PrincipalComponentsInitializer();
        TsDomain domain = dfmdata.getCurrentDomain().drop(120, 12);
        initializer.setEstimationDomain(domain);
        DynamicFactorModel model0 = initializer.initialize(dmodel, dfmdata);

        DynamicFactorModel dfm = em.initialize(model0, dfmdata);
        NumericalProcessingSpec nspec = NumericalProcessingSpec.DEFAULT_ENABLED.toBuilder()
                //                .method(NumericalProcessingSpec.Method.BFGS)
                //                .maxIntermediateIter(5)
                .build();

        DfmEstimator estimator = DfmEstimator.of(nspec, ISsfInitialization.Type.Unconditional);
        estimator.estimate(dfm, DfmEMTest.dfmdata);
        DynamicFactorModel m = estimator.getEstimatedModel();
        System.out.println(estimator.getGradient());
        
        System.out.println(m.getVar().getCoefficients());
        System.out.println(m.getVar().getInnovationsVariance());
        for (MeasurementDescriptor desc : m.getMeasurements()){
            System.out.println(desc.getVariance());
        }
            
    }

}

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
package demetra.dfm;

/**
 *
 * @author palatej
 */
@lombok.Value
@lombok.Builder(builderClassName="Builder", toBuilder=true)
public class DfmEstimationSpec implements Cloneable {
    
    @lombok.NonNull
    PrincipalComponentSpec principalComponent;
    @lombok.NonNull
    EmSpec preEm;
    @lombok.NonNull
    EmSpec postEm;
    @lombok.NonNull
    NumericalProcessingSpec processing;
    
    public static Builder builder(){
        return new Builder()
                .principalComponent(PrincipalComponentSpec.DEFAULT_ENABLED)
                .preEm(EmSpec.DEFAULT_ENABLED)
                .processing(NumericalProcessingSpec.DEFAULT_ENABLED)
                .postEm(EmSpec.DEFAULT_DISABLED);
    }
    
    public static final DfmEstimationSpec DEFAULT=builder().build();
    
    public DfmEstimationSpec disable(){
        return new Builder()
                .principalComponent(PrincipalComponentSpec.DEFAULT_DISABLED)
                .preEm(EmSpec.DEFAULT_DISABLED)
                .processing(NumericalProcessingSpec.DEFAULT_DISABLED)
                .postEm(EmSpec.DEFAULT_DISABLED)
                .build();
    }
    
    public boolean isEnabled(){
        return principalComponent.isEnabled() || preEm.isEnabled()
                ||processing.isEnabled() || postEm.isEnabled();
    }
}

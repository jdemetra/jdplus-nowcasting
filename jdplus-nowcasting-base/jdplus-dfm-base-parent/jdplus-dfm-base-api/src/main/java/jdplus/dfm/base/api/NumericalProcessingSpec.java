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
package jdplus.dfm.base.api;

/**
 *
 * @author palatej
 */
@lombok.Value
@lombok.Builder(builderClassName="Builder", toBuilder=true)
public class NumericalProcessingSpec implements Cloneable {
    
    public static enum Method{
        BFGS,
        LBFGS,
        LEVENBERGMARQUARDT
    }

    public static final int DEF_VERSION = 2, DEF_MAXITER = 1000, DEF_MAXSITER = 15,
            DEF_NITER = 5;
    public static final Boolean DEF_BLOCK = true, DEF_MIXED=true, DEF_IVAR=false;
    public static final double DEF_EPS = 1e-9;
    
    public static Builder builder(){
        return new Builder()
                 .maxIter(DEF_MAXITER)
                .maxInitialIter(DEF_MAXSITER)
                .maxIntermediateIter(DEF_NITER)
                .estimationByBlock(DEF_BLOCK)
                .mixedEstimation(DEF_MIXED)
                .independentShocks(DEF_IVAR)
                .precision(DEF_EPS)
                .method(Method.LEVENBERGMARQUARDT);
    }
    
    private boolean enabled;
    private int maxIter, maxInitialIter, maxIntermediateIter;
    private boolean estimationByBlock, mixedEstimation, independentShocks;
    private double precision;
    private Method method;

    public static final NumericalProcessingSpec DEFAULT_DISABLED=builder().build(),
            DEFAULT_ENABLED=builder().enabled(true).build();
}

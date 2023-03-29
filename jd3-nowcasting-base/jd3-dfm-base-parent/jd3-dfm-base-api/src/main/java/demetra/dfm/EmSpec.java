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
@lombok.Builder(builderClassName = "Builder", toBuilder = true)
public class EmSpec implements Cloneable {

    public static final int DEF_VERSION = 2, DEF_MAXITER = 100, DEF_MAXNUMITER = 50;
    public static final double DEF_PRECISION = 1e-9;

    private boolean enabled;
    private int version;
    private int maxIter;
    private int maxNumIter;
    private double precision;
    
    public static Builder builder(){
        return new Builder()
                .maxIter(DEF_MAXITER)
                .maxNumIter(DEF_MAXNUMITER)
                .version(DEF_VERSION)
                .precision(DEF_PRECISION);
    }
    
    public static final EmSpec DEFAULT_DISABLED=builder().build(),
            DEFAULT_ENABLED=builder().enabled(true).build();
}

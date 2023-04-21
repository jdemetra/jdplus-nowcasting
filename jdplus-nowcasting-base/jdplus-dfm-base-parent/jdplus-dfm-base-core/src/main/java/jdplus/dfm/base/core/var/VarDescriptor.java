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
package jdplus.dfm.base.core.var;

import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.api.util.Validatable;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;

/**
 * Description of a VAR component
 *
 * @author Jean Palate
 */
@lombok.Value
@lombok.Builder(builderClassName = "Builder", toBuilder = true, buildMethodName = "buildWithoutValidation")
public class VarDescriptor implements Validatable<VarDescriptor> {

    public static final double AR_DEF = .6;

    
    /**
     * Number of lags
     */
    int nlags;
    /**
     * Number of factors
     */
    int nfactors;
    /**
     * Parameters of the VAR equations The row i contains the coefficients
     * relative to the factor(i) (n=number of factors - 1) row(i): fi(t)=
     * c(i,0)f0(t-1)+...+c(i,n)fn(t-1)+... +
     * c(i,(n+1)*(k-1))f0(t-k)...+c(i,(n+1)*(k-1)+n)fn(t-k)+...
     */
    @lombok.NonNull
    Matrix coefficients;
    /**
     * Covariance matrix of the innovations
     */
    @lombok.NonNull
    private final Matrix innovationsVariance;

    public static FastMatrix defaultInnovationsCovariance(int nfactors) {
        return FastMatrix.identity(nfactors);
    }

    public static FastMatrix defaultCoefficients(int nfactors, int nlags) {
        FastMatrix c = FastMatrix.make(nfactors, nfactors * nlags);
        c.diagonal().set(AR_DEF);
        return c;
    }

    @Override
    public VarDescriptor validate() throws IllegalArgumentException {
        if (nfactors <= 0 || nlags <= 0) {
            throw new IllegalArgumentException();
        }
        if (!innovationsVariance.isSquare()
                || innovationsVariance.getRowsCount() != nfactors) {
            throw new IllegalArgumentException("Invalid innovations covariance");
        }
        if (coefficients.getRowsCount() != nfactors
                || coefficients.getColumnsCount() != nfactors * nlags) {
            throw new IllegalArgumentException("Invalid coefficents matrix");
        }
        return this;
    }

    public VarDescriptor rescaleVariance(double c) {
        if (c == 1) {
            return this;
        }
        FastMatrix var = FastMatrix.of(innovationsVariance);
        var.mul(c);
        return toBuilder()
                .innovationsVariance(var)
                .buildWithoutValidation();
    }

    public static class Coefficients {
        
        static public Coefficients of (VarDescriptor desc){
            return new Coefficients(desc);
        }

        public Coefficients(int nfactors, int nlags) {
            C = FastMatrix.make(nfactors, nlags * nfactors);
        }
        
        private Coefficients(VarDescriptor desc){
            C=FastMatrix.of(desc.coefficients);
        }

        private final FastMatrix C;

        /**
         * Gets the matrix of the var parameters corresponding to a given lag
         *
         * @param lag The lag in the var equation. Should belong to [1, nlags]
         * @return The corresponding square sub-matrix is returned. That
         * sub-matrix is a view of the underlying parameters
         */
        public FastMatrix A(int lag) {
            int nf = C.getRowsCount();
            int c0 = (lag - 1) * nf;
            return C.extract(0, nf, c0, nf);
        }

        public FastMatrix all() {
            return C;
        }

        public void setDefault() {
            C.set(0);
            C.diagonal().set(AR_DEF);
        }
    }
}

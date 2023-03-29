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
package demetra.var;

import demetra.dfm.DfmException;
import demetra.util.Validatable;
import jdplus.math.matrices.FastMatrix;

/**
 * Description of a VAR component
 * It should be noted that the structure of the descriptor can't be modified, but that the cells of the matrices
 * can be changed.
 * Make copies if need be.
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
    FastMatrix coefficients;
    /**
     * Covariance matrix of the innovations
     */
    @lombok.NonNull
    private final FastMatrix innovationsCovariance;

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
        if (nfactors <= 0 || nlags <= 0)
            throw new IllegalArgumentException();
        if (!innovationsCovariance.isSquare()
                || innovationsCovariance.getRowsCount() != nfactors) {
            throw new IllegalArgumentException("Invalid innovations covariance");
        }
        if ( coefficients.getRowsCount() != nfactors
                || coefficients.getColumnsCount() != nfactors*nlags) {
            throw new IllegalArgumentException("Invalid coefficents matrix");
        }
        return this;
    }

    public void rescaleVariance(double c) {
        innovationsCovariance.mul(c);
    }

    /**
     * Gets the matrix of the var parameters corresponding to a given lag
     *
     * @param lag The lag in the var equation. Should belong to [1, nlags]
     * @return The corresponding square sub-matrix is returned. That sub-matrix
     * is a view of the underlying parameters
     */
    public FastMatrix getA(int lag) {
        int c0 = (lag - 1) * nfactors;
        return coefficients.extract(0, nfactors, c0, nfactors);
    }

    /**
     * Sets the matrix of the var parameters corresponding to a given lag
     *
     * @param lag The lag in the var equation. Should belong to [1, nlags]
     * @param a The matrix
     */
    public void setA(int lag, FastMatrix a) {
        int n = coefficients.getRowsCount();
        for (int i = 0, j = lag - 1; i < n; ++i, j += nlags) {
            coefficients.column(j).copy(a.column(i));
        }
    }

    public void copy(VarDescriptor vdesc) {
        if (this.nfactors != vdesc.nfactors || this.nlags != vdesc.nlags) {
            throw new DfmException(DfmException.INCOMPATIBLE_DATA);
        }
        innovationsCovariance.copy(vdesc.innovationsCovariance);
        coefficients.copy(vdesc.coefficients);
    }
    
    public void setDefault(){
        coefficients.set(0);
        coefficients.diagonal().set(AR_DEF);
        innovationsCovariance.set(0);
        innovationsCovariance.diagonal().set(1);
    }

}

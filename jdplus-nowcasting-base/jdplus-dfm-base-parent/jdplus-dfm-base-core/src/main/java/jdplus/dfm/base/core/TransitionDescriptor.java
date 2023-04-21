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

import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.api.util.Validatable;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.dfm.base.core.var.VarDescriptor;

/**
 * Description of the VAR part of the model It should be noted that the order of
 * the coefficients differs from the VarDescriptor implementation
 *
 * @author palatej
 */
@lombok.Value
@lombok.Builder(builderClassName = "Builder", toBuilder = true, buildMethodName = "buildWithoutValidation")
public class TransitionDescriptor implements Validatable<TransitionDescriptor> {

    public static final double AR_DEF = .6;

    /**
     * Number of factors
     */
    int nfactors;
    /**
     * Number of lags
     */
    int nlags;
    /**
     * Parameters of the VAR equations The row i contains the coefficients
     * c(i,k) of fi(t): fi(t)= c(i,0)f0(t-1)+...+c(i,nlags-1)f0(t-nlags)+...
     * +c(i,k)fn(t-1)...+c(i,l)fn(t-nlags)). The coefficients are ordered to fit
     * the state space form (
     */
    @lombok.NonNull
    Matrix coefficients;
    /**
     * Covariance matrix of the innovations
     */
    @lombok.NonNull
    private final Matrix innovationsVariance;

    public TransitionDescriptor withDefault() {

        FastMatrix V = FastMatrix.identity(nfactors);
        FastMatrix C = FastMatrix.make(nfactors, nfactors * nlags);
        for (int i = 0; i < nfactors; ++i) {
            C.set(i, i * nlags, AR_DEF);
        }
        return builder()
                .nfactors(nfactors)
                .nlags(nlags)
                .coefficients(C)
                .innovationsVariance(V)
                .buildWithoutValidation();
    }

    public static TransitionDescriptor of(VarDescriptor v) {
        int nf = v.getNfactors(), nl = v.getNlags();
        Matrix M = v.getCoefficients();
        FastMatrix C = FastMatrix.make(nf, nf * nl);
        // M columns arrange by lags, C columns arranged by factors
        for (int l = 0, k=0; l < nl; ++l) {
            for (int f=0; f<nf; ++f, ++k){
                C.column(l+f*nl).copy(M.column(k));
            }
        }
        return builder()
                .nfactors(nf)
                .nlags(nl)
                .coefficients(C)
                .innovationsVariance(v.getInnovationsVariance())
                .buildWithoutValidation();
    }

    public VarDescriptor toVarDescriptor() {
        FastMatrix M = FastMatrix.make(nfactors, nfactors * nlags);
        // M columns arrange by lags, coefficients columns arranged by factors
        for (int f = 0, k=0; f < nfactors; ++f) {
            for (int l=0; l<nlags; ++l, ++k){
                M.column(f+l*nfactors).copy(coefficients.column(k));
            }
        }
        return VarDescriptor.builder()
                .nfactors(nfactors)
                .nlags(nlags)
                .coefficients(M)
                .innovationsVariance(innovationsVariance)
                .buildWithoutValidation();
    }
    /**
     * Gets the matrix of the var parameters corresponding to a given lag
     *
     * @param lag The lag in the var equation. Should belong to [1, nlags]
     * @return A new matrix is returned
     */
    public FastMatrix getA(int lag) {
        int n = coefficients.getRowsCount();
        FastMatrix a = FastMatrix.square(n);
        for (int i = 0, j = lag - 1; i < n; ++i, j += nlags) {
            a.column(i).copy(coefficients.column(j));
        }
        return a;
    }
    

    @Override
    public TransitionDescriptor validate() throws IllegalArgumentException {
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

    public TransitionDescriptor rescaleVariance(double c) {
        if (c == 1) {
            return this;
        }
        FastMatrix var = FastMatrix.of(innovationsVariance);
        var.mul(c);
        return toBuilder()
                .innovationsVariance(var)
                .buildWithoutValidation();
    }


}

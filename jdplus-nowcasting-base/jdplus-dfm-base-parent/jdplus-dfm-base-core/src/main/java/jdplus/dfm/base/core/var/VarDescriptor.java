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
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.math.matrices.LowerTriangularMatrix;
import jdplus.toolkit.base.core.math.matrices.SymmetricMatrix;
import org.checkerframework.checker.nullness.qual.NonNull;

/**
 * Description of a VAR component
 *
 * @author Jean Palate
 */
public class VarDescriptor {

    public static final double AR_DEF = .6;

    /**
     * Parameters of the VAR equations The row i contains the coefficients
     * relative to the factor(i) (n=number of factors - 1) row(i): fi(t)=
     * c(i,0)f0(t-1)+...+c(i,n)fn(t-1)+... +
     * c(i,(n+1)*(k-1))f0(t-k)...+c(i,(n+1)*(k-1)+n)fn(t-k)+...
     */
    private final Matrix coefficients;
    /**
     * Covariance matrix of the innovations
     */
    private final Matrix innovationsVariance;

    public VarDescriptor(@NonNull Matrix coefficients, @NonNull Matrix innovationsVariance) {

        if (coefficients.getRowsCount() != innovationsVariance.getRowsCount()
                || !innovationsVariance.isSquare()
                || coefficients.getColumnsCount() % innovationsVariance.getRowsCount() != 0) {
            throw new IllegalArgumentException("Invalid var descriptor");
        }

        this.coefficients = coefficients;
        this.innovationsVariance = innovationsVariance;
    }

    public VarDescriptor(@NonNull Matrix coefficients) {
        this.coefficients = coefficients;
        this.innovationsVariance = defaultInnovationsVariance(coefficients.getRowsCount());
    }

    /**
     * Number of lags
     *
     * @return
     */
    public int getNlags() {
        return coefficients.getColumnsCount() / coefficients.getRowsCount();
    }

    /**
     * Number of factors
     *
     * @return
     */
    public int getNfactors() {
        return coefficients.getRowsCount();
    }

    public Matrix getInnovationsVariance() {
        return innovationsVariance;
    }

    public Matrix getCoefficients() {
        return coefficients;
    }

    public static FastMatrix defaultInnovationsVariance(int nfactors) {
        return FastMatrix.identity(nfactors);
    }

    public static FastMatrix defaultCoefficients(int nfactors, int nlags) {
        FastMatrix c = FastMatrix.make(nfactors, nfactors * nlags);
        c.diagonal().set(AR_DEF);
        return c;
    }

    public VarDescriptor withDefault() {

        FastMatrix V = defaultInnovationsVariance(getNfactors());
        FastMatrix C = defaultCoefficients(getNfactors(), getNlags());
        return new VarDescriptor(V, C);
    }

    /**
     * Rescale the variance (coefficients unchanged)
     *
     * @param c
     * @return
     */
    public VarDescriptor rescaleVariance(double c) {
        if (c == 1) {
            return this;
        }
        FastMatrix var = FastMatrix.of(innovationsVariance);
        var.mul(c);
        return new VarDescriptor(coefficients, var);
    }

    /**
     * Multiply each component of the VAR by the weight given as parameter. The
     * coefficients and the innovation variance are modified accordingly
     *
     * @param w
     * @return
     */
    public VarDescriptor multiply(double[] w) {
        // covar
        int nf = getNfactors(), nc = coefficients.getColumnsCount();
        FastMatrix V = FastMatrix.of(innovationsVariance);
        FastMatrix C = FastMatrix.of(coefficients);
        for (int i = 0; i < nf; ++i) {
            V.set(i, i, 1);
            for (int j = 0; j < i; ++j) {
                V.mul(i, j, w[i] * w[j]);
            }
            for (int j = 0; j < nf; ++j) {
                if (i != j) {
                    double q = w[i] / w[j];
                    for (int c = j; c < nc; c += nf) {
                        V.mul(i, c, q);
                    }
                }
            }
        }
        SymmetricMatrix.fromLower(V);
        return new VarDescriptor(C, V);
    }

    public VarDescriptor divide(double[] w) {
        double[] iw = w.clone();
        for (int i = 0; i < w.length; ++i) {
            iw[i] = 1 / w[i];
        }
        return multiply(iw);
    }

    public static class Coefficients {

        static public Coefficients of(VarDescriptor desc) {
            return new Coefficients(desc);
        }

        public Coefficients(int nfactors, int nlags) {
            C = FastMatrix.make(nfactors, nlags * nfactors);
        }

        private Coefficients(VarDescriptor desc) {
            C = FastMatrix.of(desc.coefficients);
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
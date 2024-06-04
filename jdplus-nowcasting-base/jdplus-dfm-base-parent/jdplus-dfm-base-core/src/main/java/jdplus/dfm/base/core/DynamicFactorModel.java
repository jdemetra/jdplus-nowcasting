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

import internal.jdplus.dfm.base.core.SsfDfm;
import java.util.Collections;
import java.util.List;
import jdplus.dfm.base.api.DfmException;
import jdplus.dfm.base.core.var.VarDescriptor;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.information.GenericExplorable;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.math.matrices.LowerTriangularMatrix;
import jdplus.toolkit.base.core.math.matrices.MatrixException;
import jdplus.toolkit.base.core.math.matrices.SymmetricMatrix;
import jdplus.toolkit.base.core.ssf.multivariate.IMultivariateSsf;

/**
 *
 * @author Jean Palate
 */
@lombok.Value
@lombok.Builder(builderClassName = "Builder")
public class DynamicFactorModel implements GenericExplorable {

    /**
     * Description of the hidden VAR model
     */
    @lombok.With
    @lombok.NonNull
    private VarDescriptor var;

    /**
     * Measurement equation
     */
    @lombok.Singular
    private List<MeasurementDescriptor> measurements;

    public DynamicFactorModel rescaleVariances(double cvar) {
        if (cvar == 1) {
            return this;
        }

        Builder builder = builder().var(var.rescaleVariance(cvar));
        for (MeasurementDescriptor m : measurements) {
            builder.measurement(m.rescaleVariance(cvar));
        }
        return builder.build();
    }

    public boolean isValid() {
        for (MeasurementDescriptor mdesc : this.measurements) {
            if (mdesc.getVariance() < 0) {
                return false;
            }
        }
        FastMatrix v = FastMatrix.of(this.var.getInnovationsVariance());
        try {
            SymmetricMatrix.lcholesky(v);
        } catch (MatrixException err) {
            return false;
        }
        DfmMapping mapping = new DfmMapping(this);
        return mapping.checkBoundaries(mapping.getDefaultParameters());
    }

    /**
     * Rescale the model so that the variances of the transition shocks are
     * equal to 1. The method divides each factor by the standard deviation of
     * the corresponding transition shock and updates the different coefficients
     * accordingly.
     *
     * @return
     */
    public DynamicFactorModel normalize() {
        // scaling factors
        int nf = var.getNfactors();
        double[] w = new double[nf];
        var.getInnovationsVariance().diagonal().sqrt().copyTo(w, 0);

        Builder builder = builder().var(var.divide(w));

        // loadings
        for (MeasurementDescriptor desc : measurements) {
            double[] coefficient = desc.getCoefficient().toArray();
            for (int i = 0; i < nf; ++i) {
                if (Double.isFinite(coefficient[i])) {
                    coefficient[i] *= w[i];
                }
            }
            MeasurementDescriptor ndesc = desc.toBuilder()
                    .coefficient(DoubleSeq.of(coefficient))
                    .build();
            builder.measurement(ndesc);
        }
        return builder.build();
    }

    /**
     * Rescale the model so that the variance of the transition is I.The method
     * pre-multiplies the factor by the inverse of the Cholesky factor of the
     * covariance matrix of the transition innovations. The different
     * coefficients are updated accordingly
     *
     * @return
     * @throws DfmException is thrown when the loadings are not compatible with
     * the triangular transformation implied by Cholesky
     */
    public DynamicFactorModel lnormalize() {
        int nf = var.getNfactors(), nl = var.getNlags();
        FastMatrix V = FastMatrix.of(var.getInnovationsVariance());
        if (V.isIdentity()) {
            return this;
        }

        FastMatrix L = V.deepClone();
        SymmetricMatrix.lcholesky(L);

        // transform the var
        // f(t) = A f(t-1) + u(t)
        //L^-1*f(t) = L^-1*A*L*L^-1* f(t-1) + e(t)
        // C=L^-1*A*L <-> LC=AL
        VarDescriptor.Coefficients C = VarDescriptor.Coefficients.of(var);
        for (int i = 1; i <= nl; ++i) {
            FastMatrix A = C.A(i);
            // AL
            LowerTriangularMatrix.ML(L, A);
            // LC = (AL)
            LowerTriangularMatrix.solveLX(L, A);
        }

        Builder builder = builder().var(new VarDescriptor(C.all(), var.getInitialization()));

        // L contains the Cholesky factor
        // transform the loadings
        // y = C*f + e <-> y = (C*L)*L^-1*f+e
        // B = C*L
        // loadings
        for (MeasurementDescriptor desc : measurements) {
            double[] c = desc.getCoefficient().toArray();
            for (int i = 0; i < nf; ++i) {
                double z = 0;
                boolean nd = false;
                for (int j = i; j < nf; ++j) {
                    if (!Double.isNaN(c[j])) {
                        if (nd) {
                            throw new DfmException("Unsupported model");
                        }
                        z += c[j] * L.get(j, i);
                    } else {
                        nd = true;
                    }
                }
                if (Double.isFinite(c[i])) {
                    c[i] = z;
                }
            }
            MeasurementDescriptor ndesc = desc.toBuilder().coefficient(DoubleSeq.of(c)).build();
            builder.measurement(ndesc);
        }
        return builder.build();
    }

    public int getNfactors() {
        return var.getNfactors();
    }

    public int getNlags() {
        return var.getNlags();
    }

    public VarDescriptor getVar() {
        return var;
    }

    public List<MeasurementDescriptor> getMeasurements() {
        return Collections.unmodifiableList(measurements);
    }

    /**
     * @return Normal number of factor items in the ssf representation. The
     * state should contain: f0(t)...f0(t-n+1) f1(t)...f1(t-n+1)... where
     * n=defaultSsfBlockLength()
     */
    public int defaultSsfBlockLength() {
        int n = 1 + var.getNlags();
        for (MeasurementDescriptor desc : measurements) {
            int l = desc.getType().getLength();
            if (n < l) {
                n = l;
            }
        }
        return n;
    }

    /**
     * Number of items needed to apply the different measurements.
     * 
     * @return
     */
    public int measurementsLength() {
        int n = 0;
        for (MeasurementDescriptor desc : measurements) {
            int l = desc.getType().getLength();
            if (n < l) {
                n = l;
            }
        }
        return n;
    }

    /**
     *
     * @param nlags
     * @return
     */
    public IMultivariateSsf ssfRepresentation(int nlags) {
        return SsfDfm.of(this, nlags);
    }

    /**
     * Defines a state space with an array that contains
     * f0(t)...f0(t-blocklength+1) f1(t)...f1(t-blocklenth+1)
     * @param blockLength
     * @return 
     */
    public IMultivariateSsf ssfRepresentationWithBlockLength(int blockLength) {
        if (blockLength < defaultSsfBlockLength()) {
            throw new DfmException();
        }
        return SsfDfm.withBlockLength(this, blockLength);
    }

    /**
     *
     * @return
     */
    public int getMeasurementsCount() {
        return measurements.size();
    }

    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();
        builder.append("Loadings").append("\r\n");
        for (MeasurementDescriptor m : measurements) {
            builder.append(m).append("\r\n");
        }
        builder.append("VAR").append("\r\n");
        builder.append(var.getCoefficients());
        builder.append(var.getInnovationsVariance());
        return builder.toString();
    }

}

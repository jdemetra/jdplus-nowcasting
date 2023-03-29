/*
 * Copyright 2017 National Bank of Belgium
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

import jdplus.data.DataWindow;
import demetra.dfm.internal.SsfDfm;
import jdplus.math.matrices.FastMatrix;
import jdplus.math.matrices.LowerTriangularMatrix;
import jdplus.math.matrices.SymmetricMatrix;
import jdplus.ssf.ISsfInitialization;
import jdplus.ssf.multivariate.IMultivariateSsf;
import demetra.var.VarDescriptor;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 *
 * @author Jean Palate
 */
@lombok.Value
@lombok.Builder(builderClassName="Builder", toBuilder=true)
public class DynamicFactorModel  {

    int nlags;
    @lombok.NonNull
    VarDescriptor varDescriptor;
    ISsfInitialization.Type initialization;
    @lombok.Singular
    private List<MeasurementDescriptor> measurementDescriptors;
    private FastMatrix V0;
    
    public static Builder builder(){
        return new Builder().initialization(ISsfInitialization.Type.Unconditional);
    }

    public void rescaleVariances(double cvar) {
        for (MeasurementDescriptor m : measurementDescriptors) {
            m.rescaleVariance(cvar);
        }
        varDescriptor.rescaleVariance(cvar);
        if (V0 != null) {
            V0.mul(cvar);
        }
    }

    /**
     * Rescale the model so that the variances of the transition shocks are
     * equal to 1. The method divides each factor by the standard deviation of
     * the corresponding transition shock and updates the different coefficients
     * accordingly.
     */
    public void normalize() {
        // scaling factors
        int nl = varDescriptor.getNlags(), nf=varDescriptor.getNfactors();
        double[] w = new double[nf];
        varDescriptor.getInnovationsCovariance().diagonal().copyTo(w, 0);
        for (int i = 0; i < nf; ++i) {
            w[i] = Math.sqrt(w[i]);
        }
        if (V0 != null) {
            for (int i = 0; i < nf; ++i) {
                for (int j = 0; j < nf; ++j) {
                    V0.extract(i * nlags, nlags, j * nlags, nlags).mul(1 / (w[i] * w[j]));
                }
            }
        }
        // covar
        for (int i = 0; i < nf; ++i) {
            if (w[i] != 0) {
                varDescriptor.getInnovationsCovariance().set(i, i, 1);
                for (int j = 0; j < i; ++j) {
                    if (w[j] != 0) {
                        varDescriptor.getInnovationsCovariance().mul(i, j, 1 / (w[i] * w[j]));
                    }
                }
            }
        }
        SymmetricMatrix.fromLower(varDescriptor.getInnovationsCovariance());
        // varParams
        for (int i = 0; i < nf; ++i) {
            if (w[i] != 0) {
                DataWindow range = varDescriptor.getCoefficients().row(i).left();
                for (int j = 0; j < nf; ++j) {
                    if (w[j] != 0 && i != j) {
                        range.next(nl).mul(w[j] / w[i]);
                    }
                }
            }
        }
        // loadings
        for (MeasurementDescriptor desc : measurementDescriptors) {
            for (int i = 0; i < nf; ++i) {
                if (desc.isUsed(i)) {
                    desc.rescaleCoefficient(i, w[i]);
                }
            }
        }
    }

    /**
     * Rescale the model so that the variance of the transition is I. The method
     * pre-multiplies the factor by the inverse of the Cholesky factor of the
     * covariance matrix of the transition innovations. The different
     * coefficients are updated accordingly
     *
     * @throws A DfmException is thrown when the loadings are not compatible
     * with the triangular transformation implied by Cholesky
     */
    public void lnormalize() {
        int nf=varDescriptor.getNfactors(), nl=varDescriptor.getNlags();
        if (varDescriptor.getInnovationsCovariance().isIdentity())
            return;
        FastMatrix L = varDescriptor.getInnovationsCovariance().deepClone();
        SymmetricMatrix.lcholesky(L);
        // L contains the Cholesky factor

        // transform the loadings
        // y = C*f + e <-> y = (C*L)*L^-1*f+e
        // B = C*L
        // loadings
        for (MeasurementDescriptor desc : measurementDescriptors) {
            double[] c = desc.getCoefficients();
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
                if (desc.isUsed(i)) {
                    desc.setCoefficient(i, z);
                }
            }
        }
        // transform the var
        // f(t) = A f(t-1) + u(t)
        //L^-1*f(t) = L^-1*A*L*L^-1* f(t-1) + e(t)
        // C=L^-1*A*L <-> LC=AL
        
        for (int i = 1; i <= nl; ++i) {
            FastMatrix A = varDescriptor.getA(i);
            // AL
            LowerTriangularMatrix.ML(L, A);
            // LC = (AL)
            LowerTriangularMatrix.solveLX(L, A);
            varDescriptor.setA(i, A);
        }
        varDescriptor.getInnovationsCovariance().set(0);
        varDescriptor.getInnovationsCovariance().diagonal().set(1);
        if (V0 != null) {
            // L^-1*V*L^-1' =W <-> L(WL')=V <-> LX=V, WL'=X or LW'=X'
            for (int i = 0; i < nlags; ++i) {
                for (int j = 0; j <nlags; ++j) {
                    FastMatrix t = FastMatrix.square(nf);
                    for (int k = 0; k < nf; ++k) {
                        for (int l = 0; l < nf; ++l) {
                            t.set(k, l, this.V0.get(k * nlags + i, l * nlags + j));
                        }
                    }
                    LowerTriangularMatrix.solveLX(L, t);
                    LowerTriangularMatrix.solveLX(L, t.transpose());
                    for (int k = 0; k < nf; ++k) {
                        for (int l = 0; l < nf; ++l) {
                            V0.set(k * nlags + i, l * nlags + j, t.get(k, l));
                        }
                    }
                }
            }
        }
    }

////    /**
////     * Compacts the factors of a given models
////     *
////     * @param from The first factor to merge
////     * @param to The last factor (included) to merge
////     * @return A new model is returned. It should be re-estimated.
////     */
////    public DynamicFactorModel compactFactors(int from, int to) {
////        if (from < 0 || to < from || to >= nf_) {
////            return null;
////        }
////        if (to == from) {
////            return clone();
////        }
////        int nc = to - from;
////        DynamicFactorModel m = new DynamicFactorModel(nlags, nf_ - nc);
////        TransitionDescriptor td = new TransitionDescriptor(nf_ - nc, varDescriptor.nlags);
////        m.varDescriptor = td;
////        m.varDescriptor.covar.diagonal().set(1);
////        for (MeasurementDescriptor md : measurementDescriptors) {
////            double[] ncoeff = new double[nf_ - nc];
////            for (int i = 0; i < from; ++i) {
////                ncoeff[i] = md.coeff[i];
////            }
////            for (int i = to + 1; i < nf_; ++i) {
////                ncoeff[i - nc] = md.coeff[i];
////            }
////            boolean used = false;
////            for (int i = from; i <= to; ++i) {
////                if (!Double.isNaN(md.coeff[i])) {
////                    used = true;
////                    break;
////                }
////            }
////            if (!used) {
////                ncoeff[from] = Double.NaN;
////            }
////            m.measurementDescriptors.add(new MeasurementDescriptor(
////                    md.type, ncoeff, 1));
////        }
////        return m;
////    }
////
    /**
     * The number of lags for each factor
     *
     * @return
     */
    public int getBlockLength() {
        return nlags;
    }

    /**
     * The number of factors
     *
     * @return
     */
    public int getFactorsCount() {
        return varDescriptor.getNfactors();
    }

//    /**
//     * Changes the number of lags of each factor that is included in the model
//     *
//     * @param c The size of each block of factors (lags in [t, t-c[ belong to
//     * the model). c should larger or equal to the number of lags in the
//     * transition equation.
//     * @throws DfmException is thrown when the model is invalid (see above)
//     */
//    public void setBlockLength(int c) throws DfmException {
//        if (varDescriptor != null && c < varDescriptor.getNlags()) {
//            throw new DfmException(DfmException.INVALID_MODEL);
//        }
//        nlags = c;
//    }

    /**
     * Sets a new descriptor for the transition equation (VAR model)
     *
     * @param desc The descriptor of the transition equation
     * @return 
     * @throws DfmException is thrown when the model is invalid
     */
    public DynamicFactorModel withTransition(VarDescriptor desc) throws DfmException {
        if ( nlags < desc.getNlags()) {
            throw new DfmException(DfmException.INVALID_MODEL);
        }
        return toBuilder()
                .varDescriptor(desc)
                .build();
    }

    /**
     *
     * @return
     */
    public VarDescriptor getVarDescriptor() {
        return varDescriptor;
    }

    /**
     *
     * @return
     */
    public List<MeasurementDescriptor> getMeasurements() {
        return Collections.unmodifiableList(measurementDescriptors);
    }

    /**
     *
     * @return
     */
    public IMultivariateSsf ssfRepresentation() {
        return SsfDfm.of(varDescriptor, measurementDescriptors.toArray(MeasurementDescriptor[]::new), nlags, V0);
    }

    /**
     *
     * @return
     */
    public int getMeasurementsCount() {
        return measurementDescriptors.size();
    }

    /**
     *
     * @param init
     * @return 
     */
    public DynamicFactorModel withInitialization(ISsfInitialization.Type init) {
        Builder builder = toBuilder()
                .initialization(init);
        if (init != ISsfInitialization.Type.UserDefined) {
            builder.V0(null);
        }
        return builder.build();
    }


//    /**
//     *
//     * @return True if the model has been changed
//     */
//    public boolean validate() {
//        boolean rslt = false;
//        DfmMapping mapping = new DfmMapping(this);
//        if (!mapping.checkBoundaries(mapping.parameters())) {
//            // set default values for the VAR matrix
//            tdesc_.varParams.set(0);
//            for (int j = 0; j < nf_; ++j) {
//                tdesc_.varParams.set(j, j * tdesc_.nlags, AR_DEF);
//            }
//            rslt = true;
//        }
//
//        Matrix v = this.tdesc_.covar.clone();
//        try {
//            SymmetricMatrix.lcholesky(v);
//            return rslt;
//        } catch (MatrixException err) {
//            DataBlock d = v.diagonal().deepClone();
//            this.tdesc_.covar.set(0);
//            this.tdesc_.covar.diagonal().copy(d);
//            return true;
//        }
//
//    }
//
    public void setDefault() {
        varDescriptor.setDefault();
        for (MeasurementDescriptor m : this.measurementDescriptors) {
            m.setDefault();
        }
    }

    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();
        builder.append("Loadings").append("\r\n");
        for (MeasurementDescriptor m : measurementDescriptors) {
            builder.append(m).append("\r\n");
        }
        builder.append("VAR").append("\r\n");
        builder.append(varDescriptor.getCoefficients());
        builder.append(varDescriptor.getInnovationsCovariance());
        return builder.toString();
    }

}

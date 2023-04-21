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
package jdplus.dfm.base.core;

import internal.jdplus.dfm.base.core.SsfDfm;
import java.util.Collections;
import java.util.List;
import jdplus.dfm.base.api.DfmException;
import jdplus.dfm.base.core.var.VarDescriptor;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.core.data.DataBlock;
import jdplus.toolkit.base.core.data.DataWindow;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.math.matrices.LowerTriangularMatrix;
import jdplus.toolkit.base.core.math.matrices.MatrixException;
import jdplus.toolkit.base.core.math.matrices.SymmetricMatrix;
import jdplus.toolkit.base.core.ssf.ISsfInitialization;
import jdplus.toolkit.base.core.ssf.multivariate.IMultivariateSsf;
import lombok.NonNull;

/**
 *
 * @author Jean Palate
 */
@lombok.Value
@lombok.Builder(builderClassName = "Builder", toBuilder = true)
public class DynamicFactorModel {

    /**
     * Number of lags used in the model. Could be larger than the number of lags
     * in the var descriptor
     */
    int nlags;
    /**
     * Description of the hidden VAR model
     */
    @lombok.NonNull
    TransitionDescriptor var;
    /**
     * Initialization type
     */
    @lombok.NonNull
    ISsfInitialization.Type initialization;
    /**
     * Measurement equation
     */
    @lombok.Singular
    private List<MeasurementDescriptor> measurements;

    private DoubleSeq initialState;
    /**
     * Covariance of the initial state
     */
    private Matrix initialVariance;

    public static Builder builder() {
        return new Builder().initialization(ISsfInitialization.Type.Unconditional);
    }

    public DynamicFactorModel rescaleVariances(double cvar) {
        if (cvar == 1) {
            return this;
        }
        Builder builder = toBuilder();
        builder.clearMeasurements();
        for (MeasurementDescriptor m : measurements) {
            builder.measurement(m.rescaleVariance(cvar));
        }
        builder.var(var.rescaleVariance(cvar));
        if (initialVariance != null) {
            FastMatrix v0 = FastMatrix.of(initialVariance);
            v0.mul(cvar);
            builder.initialVariance(v0);
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
     * equal to 1.The method divides each factor by the standard deviation of
     * the corresponding transition shock and updates the different coefficients
     * accordingly.
     *
     * @return
     */
    public DynamicFactorModel normalize() {
        // scaling factors
        int nl = var.getNlags(), nf = var.getNfactors();
        double[] w = new double[nf];
        var.getInnovationsVariance().diagonal().sqrt().copyTo(w, 0);
        Builder builder = toBuilder();
        if (initialVariance != null) {
            FastMatrix V0 = FastMatrix.of(initialVariance);
            DataBlock a0 = initialState == null ? null : DataBlock.of(initialState);
            for (int i = 0; i < nf; ++i) {
                if (a0 != null)
                    a0.extract(i * nlags, nlags).mul(1/w[i]);
                for (int j = 0; j < nf; ++j) {
                    V0.extract(i * nlags, nlags, j * nlags, nlags).mul(1 / (w[i] * w[j]));
                }
            }
            builder.initialVariance(V0);
            if (a0 != null)
                builder.initialState(a0);
        }
        // covar
        TransitionDescriptor.Builder vbuilder = var.toBuilder();
        FastMatrix V = FastMatrix.of(var.getInnovationsVariance());
        for (int i = 0; i < nf; ++i) {
            if (w[i] != 0) {
                V.set(i, i, 1);
                for (int j = 0; j < i; ++j) {
                    if (w[j] != 0) {
                        V.mul(i, j, 1 / (w[i] * w[j]));
                    }
                }
            }
        }
        SymmetricMatrix.fromLower(V);
        vbuilder.innovationsVariance(V);
        
        // varParams
        FastMatrix C = FastMatrix.of(var.getCoefficients());
        for (int i = 0; i < nf; ++i) {
            if (w[i] != 0) {
                DataWindow range = C.row(i).left();
                for (int j = 0; j < nf; ++j) {
                    DataBlock b=range.next(nl);
                    if (w[j] != 0 && i != j) {
                        b.mul(w[j] / w[i]);
                    }
                }
            }
        }
        vbuilder.coefficients(C);
        builder.var(vbuilder.buildWithoutValidation());

        // loadings
        builder.clearMeasurements();
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
     * @throws DfmException is thrown when the loadings are not compatible
     * with the triangular transformation implied by Cholesky
     */
    public DynamicFactorModel lnormalize() {
        int nf = var.getNfactors(), nl = var.getNlags();
        FastMatrix V = FastMatrix.of(var.getInnovationsVariance());
        if (V.isIdentity()) {
            return this;
        }
        Builder builder = this.toBuilder();
        FastMatrix L = V.deepClone();
        SymmetricMatrix.lcholesky(L);
        // L contains the Cholesky factor

        // transform the loadings
        // y = C*f + e <-> y = (C*L)*L^-1*f+e
        // B = C*L
        // loadings
        builder.clearMeasurements();
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
        // transform the var
        // f(t) = A f(t-1) + u(t)
        //L^-1*f(t) = L^-1*A*L*L^-1* f(t-1) + e(t)
        // C=L^-1*A*L <-> LC=AL

        VarDescriptor.Coefficients C = VarDescriptor.Coefficients.of(var.toVarDescriptor());
        for (int i = 1; i <= nl; ++i) {
            FastMatrix A = C.A(i);
            // AL
            LowerTriangularMatrix.ML(L, A);
            // LC = (AL)
            LowerTriangularMatrix.solveLX(L, A);
        }

        VarDescriptor nvar=VarDescriptor.builder()
                .nfactors(var.getNfactors())
                .nlags(var.getNlags())
                .innovationsVariance(FastMatrix.identity(var.getNfactors()))
                .coefficients(C.all())
                .buildWithoutValidation();

        builder.var(TransitionDescriptor.of(nvar));
        if (initialVariance != null) {

/*          TODO
            FastMatrix V0 = FastMatrix.of(initialVariance);
            // L^-1*V*L^-1' =W <-> L(WL')=V <-> LX=V, WL'=X or LW'=X'
            for (int i = 0; i < nlags; ++i) {
                for (int j = 0; j < nlags; ++j) {
                    FastMatrix t = FastMatrix.square(nf);
                    for (int k = 0; k < nf; ++k) {
                        for (int l = 0; l < nf; ++l) {
                            t.set(k, l, this.initialVariance.get(k * nlags + i, l * nlags + j));
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
            builder.initialVariance(V);

 */
        }
        return builder.build();
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
        return var.getNfactors();
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
        if (nlags < desc.getNlags()) {
            throw new DfmException(DfmException.INVALID_MODEL);
        }
        return toBuilder()
                .var(TransitionDescriptor.of(desc))
                .build();
    }
    
    public DynamicFactorModel withBlockLength(int nlags) throws DfmException {
        if (nlags==this.nlags)
            return this;
        if (nlags < var.getNlags()) {
            throw new DfmException(DfmException.INVALID_MODEL);
        }
        return toBuilder()
                .nlags(nlags)
                .build();
    }
    

    /**
     *
     * @return
     */
    public VarDescriptor getVarDescriptor() {
        return var.toVarDescriptor();
    }

    public @NonNull TransitionDescriptor getVar() {
        return var;
    }
    /**
     *
     * @return
     */
    public List<MeasurementDescriptor> getMeasurements() {
        return Collections.unmodifiableList(measurements);
    }

    /**
     *
     * @return
     */
    public IMultivariateSsf ssfRepresentation() {
        switch (initialization){
            case Unconditional -> {
                return SsfDfm.unconditionalSsf(var, measurements.toArray(MeasurementDescriptor[]::new), nlags);
            }
             case Zero -> {
                return SsfDfm.zeroSsf(var, measurements.toArray(MeasurementDescriptor[]::new), nlags);
            }
            case UserDefined -> {
                return SsfDfm.userSsf(var, measurements.toArray(MeasurementDescriptor[]::new), initialState, initialVariance);
            }
       }
       return null;
    }

    /**
     *
     * @return
     */
    public int getMeasurementsCount() {
        return measurements.size();
    }

    /**
     *
     * @param init
     * @return
     */
    public DynamicFactorModel withInitialization(ISsfInitialization.Type init, DoubleSeq a0, Matrix V0) {
        Builder builder = toBuilder()
                .initialization(init);
        if (init != ISsfInitialization.Type.UserDefined) {
            builder.initialState(null);
            builder.initialVariance(null);
        }else{
            builder.initialVariance(V0);
            builder.initialState(a0);
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

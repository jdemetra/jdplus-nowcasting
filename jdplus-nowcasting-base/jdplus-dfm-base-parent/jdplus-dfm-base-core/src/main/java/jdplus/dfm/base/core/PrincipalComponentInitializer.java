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

import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import jdplus.dfm.base.api.PrincipalComponentSpec;
import jdplus.dfm.base.api.timeseries.TsInformationSet;
import jdplus.dfm.base.core.var.VarDescriptor;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.eco.EcoException;
import jdplus.toolkit.base.api.math.Constants;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.api.stats.AutoCovariances;
import jdplus.toolkit.base.api.timeseries.TimeSelector;
import jdplus.toolkit.base.api.timeseries.TsDomain;
import jdplus.toolkit.base.core.data.DataBlock;
import jdplus.toolkit.base.core.data.interpolation.AverageInterpolator;
import jdplus.toolkit.base.core.data.interpolation.DataInterpolator;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.math.matrices.SymmetricMatrix;
import jdplus.toolkit.base.core.math.matrices.decomposition.ISingularValueDecomposition;
import jdplus.toolkit.base.core.pca.PrincipalComponents;
import jdplus.toolkit.base.core.stats.linearmodel.LeastSquaresResults;
import jdplus.toolkit.base.core.stats.linearmodel.LinearModel;
import jdplus.toolkit.base.core.stats.linearmodel.Ols;

/**
 *
 * @author Jean Palate
 */
public class PrincipalComponentInitializer implements IDfmInitializer {

    // original figures
    private Matrix data;
    // interpolated figures
    private FastMatrix datac;
    private PrincipalComponents[] pc;
    private TsDomain domain;
    private double ns = PrincipalComponentSpec.DEF_NS;

    public TsDomain getEstimationDomain() {
        return domain;
    }

    public void setEstimationDomain(TsDomain dom) {
        domain = dom;
    }

    /**
     * Gets the minimal percentage of non missing values for determining the
     * time span of the principal components estimation.
     *
     * @return A value in ]0,1]
     */
    public double getNonMissingThreshold() {
        return ns;
    }

    public void setNonMissingThreshold(double val) {
        ns = val;
    }

    @Override
    public DynamicFactorModel initialize(DynamicFactorModel model, TsInformationSet input) {
        clear();
        // we generare the matrix corresponding to the input.
        // missing values are raughly interpolated
        if (!computeMatrix(input)) {
            return model;
        }
        // computation of the principal components on the transformed interpolated series
        // Transformations are driven by the loadings of the model

        if (!computePrincipalComponents(model)) {
            return model;
        }
        TransitionDescriptor var = computeVar(model);
        if (var == null) {
            return model;
        }
        List<MeasurementDescriptor> ndescs = computeLoadings(model);
        DynamicFactorModel nmodel = model.toBuilder()
                .clearMeasurements()
                .measurements(ndescs)
                .var(var)
                .build();
        return nmodel;
    }

    public Matrix getData() {
        return data;
    }

    public FastMatrix getInterpolatedData() {
        return datac;
    }

    public PrincipalComponents getPrincipalComponents(int block) {
        return pc[block];
    }

    private void clear() {
        data = null;
        datac = null;
        pc = null;
    }

    private boolean computeMatrix(TsInformationSet input) {
        TsDomain cdomain = this.domain;
        if (cdomain == null) {
            cdomain = searchDomain(input);
        }
        data = input.generateMatrix(cdomain);
        datac = FastMatrix.make(data.getRowsCount(), data.getColumnsCount());
        DataInterpolator interpolator = AverageInterpolator.interpolator();
        for (int i = 0; i < data.getColumnsCount(); ++i) {
            double[] col = interpolator.interpolate(data.column(i), null);
            if (col == null) {
                return false;
            }
            datac.column(i).copyFrom(col, 0);
        }
        return true;
    }

    private boolean computePrincipalComponents(DynamicFactorModel model) {
        int nb = model.getFactorsCount();
        pc = new PrincipalComponents[nb];
        for (int i = 0; i < nb; ++i) {
            FastMatrix x = prepareDataForComponent(model, i);
            pc[i] = new PrincipalComponents();
            pc[i].process(x);
        }
        return true;
    }

    /**
     * Creates the data used for the computation of the principal components
     * analysis
     *
     * @param model
     * @param cmp The considered factor
     * @return
     */
    private FastMatrix prepareDataForComponent(DynamicFactorModel model, int cmp) {
        // Keep only the concerned series
        int np = 0;
        for (MeasurementDescriptor desc : model.getMeasurements()) {
            if (!Double.isNaN(desc.getCoefficient(cmp))) {
                ++np;
            }
        }
        FastMatrix m = FastMatrix.make(datac.getRowsCount(), np);
        // Copy the series and correct them by the effect of the previous factors
        np = 0; // the position of the series in the matrix
        int s = 0; // its position in the model
        for (MeasurementDescriptor desc : model.getMeasurements()) {
            if (!Double.isNaN(desc.getCoefficient(cmp))) {
                m.column(np).copy(datac.column(s));
                for (int j = 0; j < cmp; ++j) {
                    if (!Double.isNaN(desc.getCoefficient(j))) {
                        ISingularValueDecomposition svd = pc[j].getSvd();
                        double scaling = pc[j].getScaling();
                        double l = -svd.S().get(0) * svd.V().get(searchPos(model, s, j), 0) / scaling;
                        m.column(np).addAY(l, svd.U().column(0));
                    }
                }
                ++np;
            }
            ++s;
        }
        return m;
    }

    /**
     * Searches the position of variable v in the singular value decomposition k
     *
     * @param model
     * @param v
     * @param k
     * @return
     */
    private int searchPos(DynamicFactorModel model, int v, int k) {
        int s = -1, q = 0;
        for (MeasurementDescriptor desc : model.getMeasurements()) {
            if (!Double.isNaN(desc.getCoefficient(k))) {
                ++s;
            }
            if (q == v) {
                break;
            }
            ++q;
        }
        return s;
    }

    private TransitionDescriptor computeVar(DynamicFactorModel model) {
        TransitionDescriptor var = model.getVar();
        int nl = var.getNlags(), nb = model.getFactorsCount();
        DoubleSeq[] f = new DoubleSeq[nb];
        DoubleSeq[] e = new DoubleSeq[nb];
        FastMatrix M = FastMatrix.make(data.getRowsCount() - nl, nl * nb);
        int c = 0;
        for (int i = 0; i < nb; ++i) {
            DataBlock cur = pc[i].getFactor(0);
            f[i] = cur.drop(nl, 0);
            for (int j = 1; j <= nl; ++j) {
                M.column(c++).copy(cur.drop(nl - j, j));
            }
        }
        LinearModel.Builder regmodel = LinearModel.builder();
        for (int j = 0; j < M.getColumnsCount(); ++j) {
            regmodel.addX(M.column(j));
        }
        FastMatrix C = FastMatrix.make(nb, nb*nl);
        for (int i = 0; i < nb; ++i) {
            regmodel.y(f[i]);
            try {
                LeastSquaresResults ols = Ols.compute(regmodel.build());
                C.row(i).copy(ols.getLikelihood().coefficients());
                e[i] = ols.residuals();
            } catch (EcoException ex) {
                return null;
            }
        }
        FastMatrix V = FastMatrix.of(var.getInnovationsVariance());

        for (int i = 0; i < nb; ++i) {
            for (int j = 0; j <= i; ++j) {
                V.set(i, j, AutoCovariances.covarianceWithZeroMean(e[i], e[j]));
            }
        }
        SymmetricMatrix.fromLower(V);
        return var.toBuilder()
                .coefficients(C)
                .innovationsVariance(V)
                .buildWithoutValidation();
    }

    private List<MeasurementDescriptor> computeLoadings(DynamicFactorModel model) {
        // creates the matrix of factors
        int nb = model.getFactorsCount(), blen = model.getBlockLength();
        FastMatrix M = FastMatrix.make(data.getRowsCount() - (blen - 1), nb * blen);
        for (int i = 0, c = 0; i < nb; ++i) {
            DataBlock cur = pc[i].getFactor(0);
            for (int j = 0; j < blen; ++j) {
                M.column(c++).copy(cur.drop(blen - 1 - j, j));
            }
        }
        int v = 0;
        List<MeasurementDescriptor> ndescs = new ArrayList<>();
        for (MeasurementDescriptor desc : model.getMeasurements()) {
            DataBlock y = datac.column(v++).drop(blen - 1, 0);
            if (y.isZero(Constants.getEpsilon())) {
                ndescs.add(desc.withVariance(1));
            } else {
                MeasurementDescriptor.Builder builder = desc.toBuilder();
                LinearModel.Builder regmodel = LinearModel.builder();
                regmodel.y(y);
                for (int j = 0; j < nb; ++j) {
                    if (!Double.isNaN(desc.getCoefficient(j))) {
                        double[] x = new double[y.length()];
                        int s = j * blen, l = desc.getType().getLength();
                        for (int r = 0; r < x.length; ++r) {
                            x[r] = desc.getType().dot(M.row(r).extract(s, l));
                        }
                        regmodel.addX(DoubleSeq.of(x));
                    }
                }
                try {
                    LeastSquaresResults ols = Ols.compute(regmodel.build());
                    double[] b = ols.getCoefficients().toArray();
                    double[] c = desc.getCoefficient().toArray();
                    for (int i = 0, j = 0; j < nb; ++j) {
                        if (!Double.isNaN(c[j])) {
                            c[j] = b[i++];
                        }
                    }
                    builder.coefficient(DoubleSeq.of(c))
                            .variance(ols.getErrorMeanSquares());
                } catch (EcoException ex) {
                    builder.variance(1);
                }
                ndescs.add(builder.build());
            }
        }
        return ndescs;
    }

    private TsDomain searchDomain(TsInformationSet input) {
        int n = input.getSeriesCount();
        LocalDateTime[] start = new LocalDateTime[n];
        LocalDateTime[] end = new LocalDateTime[n];
        for (int i = 0; i < n; ++i) {
            TsDomain cur = input.series(i).cleanExtremities().getDomain();
            start[i] = cur.getStartPeriod().start();
            end[i] = cur.getEndPeriod().start();
        }
        Arrays.sort(start);
        Arrays.sort(end);
        int t = (int) ((n - 1) * ns);
        TimeSelector sel = TimeSelector.between(start[t], end[n - 1 - t]);
        return input.getCurrentDomain().select(sel);
    }

}

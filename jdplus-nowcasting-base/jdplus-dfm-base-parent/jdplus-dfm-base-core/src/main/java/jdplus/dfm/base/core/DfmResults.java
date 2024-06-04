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

import jdplus.dfm.base.api.timeseries.TsInformationSet;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.information.GenericExplorable;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.core.data.DataBlock;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.ssf.StateStorage;
import jdplus.toolkit.base.core.ssf.multivariate.IMultivariateSsf;
import jdplus.toolkit.base.core.stats.likelihood.Likelihood;

/**
 *
 * @author palatej
 */
@lombok.Value
@lombok.Builder(builderClassName = "Builder")
public class DfmResults implements GenericExplorable {

    TsInformationSet dfmData;
    DynamicFactorModel dfm;
    StateStorage smoothedStates;
    DoubleSeq gradient;
    Matrix hessian;
    Likelihood likelihood;
    boolean standardizedInput;

    DoubleSeq sampleMean;
    DoubleSeq sampleStDev;
    Matrix inputData;
    Matrix transformedInputData;
    Matrix factors;
    Matrix factorsStdErr;
    Matrix residuals;
    Matrix residualsStandardized;
    Double logLikelihood;

    public Matrix forecastsT(int nf) {

        if (nf <= 0) {
            nf = dfmData.series(0).getAnnualFrequency();
        }
        int nPeriod = dfmData.getCurrentDomain().getLength();
        int nOutCalc = smoothedStates.item(0).length() - nPeriod;
        nf = Math.min(nf, nOutCalc);

        int nPeriodExt = dfmData.getCurrentDomain().getLength() + nf;
        int neq = dfm.getMeasurementsCount();

        FastMatrix fcastsT = FastMatrix.make(nPeriodExt, neq);
        IMultivariateSsf ssf = dfm.ssfRepresentation(0);

        for (int k = 0; k < neq; ++k) {
            // Position of the last observed value
            int nlk = 0;
            for (int l = nPeriod - 1; l >= 0; --l) {
                boolean isLast = Double.isFinite(dfmData.series(k).get(l).getValue());
                if (isLast) {
                    nlk = l;
                    break;
                }
            }
            // Fill observed values
            int i = 0;
            while (i <= nlk) {
                fcastsT.set(i, k, dfmData.series(k).getValue(i));
                ++i;
            }
            // Forecasts
            DataBlock Zk = DataBlock.make(ssf.getStateDim());
            ssf.measurements().loading(k).Z(k, Zk);
            int nfcsts = nPeriodExt - (nlk + 1);
            for (int h = 0; h < nfcsts; ++h) {
                DataBlock ahk = smoothedStates.a(nlk + 1 + h);
                fcastsT.set(nlk + 1 + h, k, Zk.dot(ahk));
            }
        }
        return fcastsT;
    }

    public Matrix forecasts(int nf) {
        Matrix fcastsT = forecastsT(nf);
        FastMatrix fcasts = FastMatrix.make(fcastsT.getRowsCount(), fcastsT.getColumnsCount());

        for (int j = 0; j < fcasts.getColumnsCount(); ++j) {
            DoubleSeq sj = fcastsT.column(j);
            if (!standardizedInput) {
                sj = sj.times(sampleStDev.get(j));
                sj = sj.plus(sampleMean.get(j));
            }
            fcasts.column(j).add(sj);
        }
        return (fcasts);
    }

    public Matrix forecastsTStDev(int nf) {
        if (nf <= 0) {
            nf = dfmData.series(0).getAnnualFrequency();
        }
        int nPeriod = dfmData.getCurrentDomain().getLength();
        int nOutCalc = smoothedStates.item(0).length() - nPeriod;
        nf = Math.min(nf, nOutCalc);

        int nPeriodExt = dfmData.getCurrentDomain().getLength() + nf;
        int neq = dfm.getMeasurementsCount();
        FastMatrix fcastsTStDev = FastMatrix.make(nPeriodExt, neq);
        IMultivariateSsf ssf = dfm.ssfRepresentation(0);
        
        for (int k = 0; k < neq; ++k) {
            DataBlock Zk = DataBlock.make(ssf.getStateDim());
            ssf.measurements().loading(k).Z(k, Zk);
            DoubleSeq ZPZk = smoothedStates.zvariance(Zk);
            double H = dfm.getMeasurements().get(k).getVariance();
            DoubleSeq Fk = ZPZk.plus(H);
            // in-sample
            for (int i = 0; i < nPeriod; ++i) {
                if (Double.isFinite(dfmData.series(k).get(i).getValue())) {
                    fcastsTStDev.set(i, k, 0);
                } else {
                    fcastsTStDev.set(i, k, Math.sqrt(Fk.get(i)));
                }
            }
            // out-of-sample
            for (int h = nPeriod; h < nPeriodExt; ++h) {
                fcastsTStDev.set(h, k, Math.sqrt(Fk.get(h)));
            }
        }
        return fcastsTStDev;
    }

    public Matrix forecastsStDev(int nf) {
        Matrix fcastsTStDev = forecastsTStDev(nf);
        FastMatrix fcastsStDev = FastMatrix.make(fcastsTStDev.getRowsCount(), fcastsTStDev.getColumnsCount());

        for (int i = 0; i < fcastsStDev.getRowsCount(); ++i) {
            for (int j = 0; j < fcastsStDev.getColumnsCount(); ++j) {
                double sij = fcastsTStDev.get(i, j);
                if (sij == 0) {
                    fcastsStDev.set(i, j, 0);
                } else {
                    fcastsStDev.set(i, j, (sij * sampleStDev.get(j)));
                }
            }
        }
        return (fcastsStDev);
    }

}

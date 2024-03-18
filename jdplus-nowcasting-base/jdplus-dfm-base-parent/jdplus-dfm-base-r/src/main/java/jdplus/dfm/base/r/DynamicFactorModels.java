/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package jdplus.dfm.base.r;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import jdplus.dfm.base.api.MeasurementType;
import jdplus.dfm.base.api.timeseries.TsInformationSet;
import jdplus.dfm.base.api.timeseries.TsInformationUpdates;
import jdplus.dfm.base.core.DfmEM;
import jdplus.dfm.base.core.DfmEstimates;
import jdplus.dfm.base.core.DfmEstimator;
import jdplus.dfm.base.core.DfmKernel;
import jdplus.dfm.base.core.DfmNews;
import jdplus.dfm.base.core.DfmProcessor;
import jdplus.dfm.base.core.DfmResults;
import jdplus.dfm.base.core.DfmResultsNews;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.dfm.base.core.DynamicFactorModel;
import jdplus.dfm.base.core.IDfmMeasurement;
import jdplus.dfm.base.core.MeasurementDescriptor;
import jdplus.dfm.base.core.PrincipalComponentsInitializer;
import jdplus.dfm.base.core.var.VarDescriptor;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsDomain;
import jdplus.toolkit.base.api.timeseries.TsPeriod;
import jdplus.toolkit.base.core.math.functions.levmar.LevenbergMarquardtMinimizer;
import jdplus.toolkit.base.core.math.functions.ssq.ProxyMinimizer;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.ssf.ISsfInitialization;
import jdplus.toolkit.base.r.timeseries.TsUtility;

/**
 *
 * @author LEMASSO
 */
@lombok.experimental.UtilityClass
public class DynamicFactorModels { 
    
    private TsInformationSet prepareInput(Matrix data, int freq, int[] start, boolean standardized, double[] fixedSampleMean, double[] fixedStDev, 
            DfmResults.Builder builder) {
        List<TsData> input = new ArrayList<>();
        double[] sampleMean = new double[data.getColumnsCount()];
        double[] sdDev = new double[data.getColumnsCount()];

        for (int j = 0; j < data.getColumnsCount(); ++j) {
            DoubleSeq sj = data.column(j);
            double m = fixedSampleMean == null ? sj.averageWithMissing() : fixedSampleMean[j];
            double e;
            if(fixedStDev == null){
                double e2 = sj.ssqcWithMissing(m) / sj.count(Double::isFinite);
                e = Math.sqrt(e2);
            } else{
                e = fixedStDev[j];
            }         

            if (!standardized) {
                sj = sj.plus(-m);
                sj = sj.times(1 / e);
            }

            input.add(TsUtility.of(freq, start[0], start[1], sj.toArray()));
            sampleMean[j] = m;
            sdDev[j] = e;
        }
        TsInformationSet tinput = new TsInformationSet(input);
        
        FastMatrix dataT = FastMatrix.of(data);
        for(int j = 0; j < tinput.getSeriesCount(); ++j){
            dataT.column(j).set(tinput.series(j).getValues());
        }
        if (builder != null) {
            builder.sampleMean(DoubleSeq.of(sampleMean))
                   .sampleStDev(DoubleSeq.of(sdDev))
                   .inputData(data)
                   .standardizedInput(standardized)
                   .transformedInputData(dataT);
        }        
        return tinput;
    } 
    
    public Matrix getFactors(DfmProcessor processor, int nPeriodExtended, int nFactors, int blockLength){
        FastMatrix factors = FastMatrix.make(nPeriodExtended, nFactors);       
        for (int j = 0; j < nFactors; ++j) {
            factors.column(j).add(processor.getSmoothingResults().getComponent(blockLength*j));  
        } 
        return factors;
    }
    
    public Matrix getFactorsSdErr(DfmProcessor processor, int nPeriodExtended, int nFactors, int blockLength){
        FastMatrix factorsSdErr = FastMatrix.make(nPeriodExtended, nFactors);   
        for (int j = 0; j < nFactors; ++j) {
            factorsSdErr.column(j).add(processor.getSmoothingResults().getComponentVariance(blockLength*j).sqrt());
        } 
        return factorsSdErr;
    }
    
    public Matrix getResiduals(DfmProcessor processor, TsInformationSet dfmData, int nPeriod, int nEq){
        FastMatrix vt = FastMatrix.make(nPeriod, nEq);
        vt.set(Double.NaN);
        for(int i = 0; i < nPeriod; ++i){
           int k = 0; 
           for(int j = 0; j < nEq; ++j){
               if(Double.isFinite(dfmData.series(j).get(i).getValue())){
                   vt.set(i, j, processor.getFilteringResults().get(i).getE().get(k));
                   ++k;
               }
           }    
        }
        return vt;
    }
    
    public Matrix getStandardizedResiduals(DfmProcessor processor, TsInformationSet dfmData, int nPeriod, int nEq){
        FastMatrix et = FastMatrix.make(nPeriod, nEq);
        et.set(Double.NaN);
        for(int i = 0; i < nPeriod; ++i){
           int k = 0; 
           for(int j = 0; j < nEq; ++j){
               if(Double.isFinite(dfmData.series(j).get(i).getValue())){
                   et.set(i, j, processor.getFilteringResults().get(i).getU().get(k));
                   ++k;
               }
           }    
        }
        return et;
    }
    
    public DynamicFactorModel getDfm(DfmEstimates dfmEst){
        return dfmEst.getDfm();
    }
    
    // model from scratch -> NOT USED ANYMORE IN rjd3nowcasting v2+
    public DynamicFactorModel model(int nfactors, int nlags, String[] factorType, Matrix factorLoaded, String varInit, double[] mVariance) {

        if (factorLoaded.getRowsCount() != factorType.length) {
            throw new IllegalArgumentException("The number of rows of factorLoaded should match the length of factor Type");
        }
        if (factorLoaded.getColumnsCount() != nfactors) {
            throw new IllegalArgumentException("The number of columns of factorLoaded should match nfactors");
        }

        ISsfInitialization.Type varInitType = ISsfInitialization.Type.Unconditional;
        if (varInit.equals("Zero")) {
            varInitType = ISsfInitialization.Type.Zero;
        }
        VarDescriptor varDesc = VarDescriptor.defaultVar(nfactors, nlags, varInitType);
   
        FastMatrix factorLoadedM = FastMatrix.make(factorLoaded.getRowsCount(), factorLoaded.getColumnsCount());
        List<MeasurementDescriptor> mDescs = new ArrayList<>();
        for (int i = 0; i < factorType.length; ++i) {
            IDfmMeasurement factorTransf;
            factorTransf = switch (factorType[i]) {
                case "YoY" ->
                    IDfmMeasurement.measurement(MeasurementType.YoY);
                case "Q" ->
                    IDfmMeasurement.measurement(MeasurementType.Q);
                default ->
                    IDfmMeasurement.measurement(MeasurementType.M);
            };
            for (int j = 0; j < nfactors; ++j) {
                factorLoadedM.set(i, j, factorLoaded.get(i, j) == 0 ? Double.NaN : 1); // if not used, defined as NaN in class DynamicFactorModel
            }
            MeasurementDescriptor mdesc = MeasurementDescriptor.builder()
                    .type(factorTransf)
                    .coefficient(factorLoadedM.row(i))
                    .variance(mVariance == null ? 1 : mVariance[i])
                    .build();
            mDescs.add(mdesc);
        }

        DynamicFactorModel dfm = DynamicFactorModel.builder()
                .measurements(mDescs)
                .var(varDesc)
                .build();

        return dfm;
    }
    
    // model with pre-initialized values of the parameters  -> used from rjd3nowcasting v2+
    public DynamicFactorModel model(int nfactors, int nlags, String[] factorType, Matrix factorLoaded, String varInit,
                                     Matrix varCoefficients, Matrix varInnovationsVariance, Matrix mCoefficients, double[] mIdiosyncraticErrVariance) {
        
        int nSeries = factorType.length;
        int nCoefVar = nfactors * nlags;   
        
        if (factorLoaded.getRowsCount() != nSeries) {
            throw new IllegalArgumentException("The number of rows of factorLoaded should match the length of factor Type");
        }
        if (factorLoaded.getColumnsCount() != nfactors) {
            throw new IllegalArgumentException("The number of columns of factorLoaded should match nfactors");
        } 
        if(varCoefficients != null && (varCoefficients.getRowsCount() != nfactors || varCoefficients.getColumnsCount() != nCoefVar)){
            throw new IllegalArgumentException("The dimension of the matrix related to the VAR coefficients is not consistent with the number of factors and/or lags.");
        }
        if(varInnovationsVariance != null && (varInnovationsVariance.getRowsCount() != nfactors || !varInnovationsVariance.isSquare())){
            throw new IllegalArgumentException("The dimension of the matrix related to the VAR errors variance is not consistent with the number of factors.");
        }
        if(mCoefficients != null && (mCoefficients.getRowsCount() != nSeries || mCoefficients.getColumnsCount() !=  nfactors)){
            throw new IllegalArgumentException("The dimension of the matrix related to the measurement coefficients is not consistent with the number of series (length of factorType) and/or the number of factors.");
        }
        if(mIdiosyncraticErrVariance != null && mIdiosyncraticErrVariance.length != nSeries){
            throw new IllegalArgumentException("The length of the vector related to the measurement errors variance is not consistent with the number of series (length of factorType).");
        }
        
        varCoefficients = varCoefficients == null ? VarDescriptor.defaultCoefficients(nfactors, nlags) : varCoefficients;
        varInnovationsVariance = varInnovationsVariance == null ? VarDescriptor.defaultInnovationsVariance(nfactors) : varInnovationsVariance; 
        ISsfInitialization.Type varInitType = varInit.equals("Zero") ? ISsfInitialization.Type.Zero : ISsfInitialization.Type.Unconditional;
        VarDescriptor varDesc = new VarDescriptor(varCoefficients, varInnovationsVariance, varInitType);
        
        if(mCoefficients == null){
            FastMatrix tmp = FastMatrix.make(factorLoaded.getRowsCount(), factorLoaded.getColumnsCount());
            tmp.set(1);
            mCoefficients = tmp;
        }
        if(mIdiosyncraticErrVariance == null){
           mIdiosyncraticErrVariance = new double[nSeries]; 
           Arrays.fill(mIdiosyncraticErrVariance, 1);
        } 
        FastMatrix mCoefficientsFormatted = FastMatrix.make(factorLoaded.getRowsCount(), factorLoaded.getColumnsCount());
        List<MeasurementDescriptor> mDescs = new ArrayList<>();
        for (int i = 0; i < nSeries; ++i) {
            IDfmMeasurement factorTransf;
            factorTransf = switch (factorType[i]) {
                case "YoY" ->
                    IDfmMeasurement.measurement(MeasurementType.YoY);
                case "Q" ->
                    IDfmMeasurement.measurement(MeasurementType.Q);
                default ->
                    IDfmMeasurement.measurement(MeasurementType.M);
            };
            for (int j = 0; j < nfactors; ++j) {
                mCoefficientsFormatted.set(i, j, factorLoaded.get(i, j) == 0 ? Double.NaN : mCoefficients.get(i, j));
            }
            MeasurementDescriptor mdesc = MeasurementDescriptor.builder()
                    .type(factorTransf)
                    .coefficient(mCoefficientsFormatted.row(i))
                    .variance(mIdiosyncraticErrVariance[i])
                    .build();
            mDescs.add(mdesc);
        }

        DynamicFactorModel dfm = DynamicFactorModel.builder()
                .measurements(mDescs)
                .var(varDesc)
                .build();
        if (!dfm.isValid()) {
            throw new IllegalArgumentException("Invalid DFM");
        } 
        
        return dfm;
    }
    
    public DfmEstimates estimate_PCA(DynamicFactorModel dfmModel, Matrix data, int freq, int[] start, boolean standardized){
                
        TsInformationSet dfmData = prepareInput(data, freq, start, standardized, null, null, null);
        
        PrincipalComponentsInitializer initializer = new PrincipalComponentsInitializer();
        DynamicFactorModel dfm = initializer.initialize(dfmModel, dfmData);
        
        DfmEstimates dfmE = DfmEstimates.builder()
                .dfm(dfm)
                .build();
        
        return dfmE; 
    }
      
    public DfmEstimates estimate_EM(DynamicFactorModel dfmModel, Matrix data, int freq, int[] start, boolean standardized, 
            boolean pcaInit, int maxIter, double eps) {
        
        TsInformationSet dfmData = prepareInput(data, freq, start, standardized, null, null, null);
        
        DynamicFactorModel dfm0;
        if (pcaInit) {
            PrincipalComponentsInitializer initializer = new PrincipalComponentsInitializer();
            DynamicFactorModel pcaModel = initializer.initialize(dfmModel, dfmData);
            if (pcaModel.isValid()) {
                dfm0 = pcaModel;
            }else{
                dfm0 = dfmModel;
            }
        }else{
            dfm0 = dfmModel;
        }
        
        DfmEM em = DfmEM.builder()
                .maxIter(maxIter)
                .precision(eps)
                .build();   
        DynamicFactorModel dfm = em.initialize(dfm0, dfmData);
        
        DfmEstimates dfmE = DfmEstimates.builder()
                .dfm(dfm)
                .ll(em.getFinalLogLikelihood())
                .build();
        
        return dfmE; 
    }

    public DfmEstimates estimate_ML(DynamicFactorModel dfmModel, Matrix data, int freq, int[] start, boolean standardized,
            boolean pcaInit, boolean emInit, int emInitmaxIter, double emInitEps, int maxIter, int maxBlockIter, 
            int simplModelIter, boolean independentVARShocks, boolean mixedEstimation, double eps) {
            
        TsInformationSet dfmData = prepareInput(data, freq, start, standardized, null, null, null);
         
        DynamicFactorModel dfm0 = dfmModel;
        if (pcaInit) {
            PrincipalComponentsInitializer pcaInitializer = new PrincipalComponentsInitializer();
            dfm0 = pcaInitializer.initialize(dfmModel, dfmData);
        }
        if (emInit) { 
            DfmEM emInitializer = DfmEM.builder()
                .maxIter(emInitmaxIter)
                .precision(emInitEps)
                .build(); 
            dfm0 = emInitializer.initialize(dfm0, dfmData);
        }    
      
        DfmEstimator estimator = DfmEstimator.builder()
                .minimizer(ProxyMinimizer.builder(
                        LevenbergMarquardtMinimizer.builder()))
                .maxInitialIter(simplModelIter)
                .maxBlockIterations(maxBlockIter)
                .maxIterations(maxIter)
                .independentVarShocks(independentVARShocks)
                .mixed(mixedEstimation)
                .precision(eps)
                .build();
        

        estimator.estimate(dfm0, dfmData);
        
        DfmEstimates dfmE = DfmEstimates.builder()
                .dfm(estimator.getEstimatedModel())
                .ll(estimator.getLikelihood().logLikelihood())
                .hessian(estimator.getHessian())
                .gradient(estimator.getGradient().toArray())
                .hasConverged(estimator.hasConverged())
                .build();
        
        return dfmE;
    }
    
    public DfmResults process(DynamicFactorModel dfmModel, Matrix data, int freq, int[] start, boolean standardized, int nForecasts){
        
        DfmResults.Builder builder = DfmResults.builder();
        
        TsInformationSet dfmData = prepareInput(data, freq, start, standardized, null, null, builder);
        
        DfmProcessor processor = DfmProcessor.builder().build();
        TsPeriod fLast = dfmData.getCurrentDomain().getEndPeriod().plus(nForecasts);
        TsInformationSet dfmDataExtended = dfmData.extendTo(fLast.start().toLocalDate());
        processor.process(dfmModel, dfmDataExtended);
        
        int nPeriod = dfmData.getCurrentDomain().getLength();
        int nPeriodExt = dfmDataExtended.getCurrentDomain().getLength();
        int nf = dfmModel.getNfactors();
        int neq =dfmModel.getMeasurementsCount();
        int blockLength = dfmModel.defaultSsfBlockLength();
        Matrix factors = getFactors(processor, nPeriodExt, nf, blockLength);
        Matrix factorsSdErr = getFactorsSdErr(processor, nPeriodExt, nf, blockLength);
        Matrix residuals = getResiduals(processor, dfmData, nPeriod, neq);
        Matrix residualsStandardized = getStandardizedResiduals(processor, dfmData, nPeriod, neq);
         
        return builder
                .dfmData(dfmData)
                .dfm(dfmModel)
                .smoothedStates(processor.getSmoothingResults())
                .factors(factors)
                .factorsStdErr(factorsSdErr)
                .residuals(residuals)
                .residualsStandardized(residualsStandardized)
                .build();    
    }
        
    public DfmResultsNews computeNews(int targetSeries, DynamicFactorModel dfm, Matrix oldData, Matrix newData, int freq, int[] start, boolean standardized, int nForecasts){ 
        
        if(oldData.getColumnsCount() != newData.getColumnsCount()){
            throw new IllegalArgumentException("The number of columns in old and new dataset is different!");
        }
        
        // Preprocessing: same transformation must be applied to oldData and newData   
        TsInformationSet dfmOldData = prepareInput(oldData, freq, start, standardized, null, null, null);  
        TsInformationSet dfmNewData;
        int nSeries = oldData.getColumnsCount();
        double[] oldDataSampleMean = new double[nSeries];
        double[] oldDataSdDev = new double[nSeries];
        
        if(!standardized){
            for (int j = 0; j < nSeries; ++j) {
                DoubleSeq sj = oldData.column(j);
                double m = sj.averageWithMissing();
                double e2 = sj.ssqcWithMissing(m) / sj.count(Double::isFinite);
                double e = Math.sqrt(e2);
                oldDataSampleMean[j] = m;
                oldDataSdDev[j] = e;
            }   
            dfmNewData = prepareInput(newData, freq, start, standardized, oldDataSampleMean, oldDataSdDev, null);
        } else{
            dfmNewData = prepareInput(newData, freq, start, true, null, null, null);
        }
        
        TsPeriod fcstLast = dfmNewData.getCurrentDomain().getEndPeriod().plus(nForecasts);
        TsInformationSet dfmOldDataExtended = dfmOldData.extendTo(fcstLast.start().toLocalDate());
        TsInformationSet dfmNewDataExtended = dfmNewData.extendTo(fcstLast.start().toLocalDate());
        
        // Computation of news (pm news = newSet - revisedSet)
        DfmNews newsData = new DfmNews(dfm);
        boolean isUpdates = newsData.process(dfmOldDataExtended, dfmNewDataExtended);
        if(!isUpdates){
            throw new IllegalArgumentException("No updates between the two datasets!");
        }
        
        // Output (for a given series of interest) 
        
        // News
        int nn = newsData.news().length(); 
        int[] index = new int[nn];
        String[] per = new String[nn];
        double[] expT, obsT, newsT, exp, obs, news;
        expT = new double[nn];
        obsT = new double[nn];
        newsT = new double[nn];
        exp = new double[nn];
        obs = new double[nn];
        news = new double[nn];
          
        List<TsInformationUpdates.Update> updates = newsData.newsDetails().news();        
        for (int k = 0; k < nn; ++k) {      
            index[k] = updates.get(k).getSeries();
            TsPeriod pi = updates.get(k).getPeriod();
            per[k] = pi.display();
            expT[k] = updates.get(k).getForecast();
            obsT[k] = updates.get(k).getObservation();
            newsT[k] = updates.get(k).getNews();
            
            // Revert transformation if any
            if(!standardized){
                double meank = oldDataSampleMean[index[k]];
                double sdk = oldDataSdDev[index[k]];    
                exp[k] = (expT[k] * sdk) + meank;
                obs[k] = (obsT[k] * sdk) + meank;
                news[k] = obs[k] - exp[k];        
            }else{
                exp[k] = expT[k];
                obs[k] = obsT[k]; 
                news[k] = newsT[k]; 
            }
        }
        
        // Weights & impacts
        int nf =  nForecasts;
        FastMatrix sWeightsT, sImpactsT, sWeights, sImpacts;
        sWeightsT = FastMatrix.make(nn, nf); 
        sImpactsT = FastMatrix.make(nn, nf);
        sWeights = FastMatrix.make(nn, nf); 
        sImpacts = FastMatrix.make(nn, nf);
        TsPeriod endPerTarget = dfmNewData.series(targetSeries).cleanExtremities().getEnd();
        double meanTarget = oldDataSampleMean[targetSeries];
        double sdTarget = oldDataSdDev[targetSeries];
        
        for (int j = 0; j < nf; ++j){
            DoubleSeq bTj = newsData.weights(targetSeries, endPerTarget.plus(j));
            sWeightsT.column(j).set(bTj);
            
            for (int i = 0; i < nn; ++i){
                sImpactsT.set(i, j, bTj.get(i) * newsT[i]);
            }
        }
        
        if(!standardized){
            for (int i = 0; i < nn; ++i){
                DoubleSeq bTi = DoubleSeq.of(sWeightsT.row(i).toArray());
                sWeights.row(i).set(bTi.times(sdTarget/oldDataSdDev[index[i]]));
                DoubleSeq bi = DoubleSeq.of(sWeights.row(i).toArray());
                sImpacts.row(i).set(bi.times(news[i]));
            }
        }else{
            sWeights = sWeightsT;
            sImpacts = sImpactsT;    
        } 
        
        // Forecasts     
        String[] fcstPeriods = new String[nf]; 
        double[] oldFcstsT, oldFcsts, revFcstsT, revFcsts, newFcstsT, newFcsts;
        oldFcstsT = new double[nf];
        revFcstsT = new double[nf];
        newFcstsT = new double[nf];
        oldFcsts = new double[nf];
        revFcsts = new double[nf];
        newFcsts = new double[nf];
              
        for (int j = 0; j < nf; ++j){
            fcstPeriods[j] = endPerTarget.plus(j).display();
            oldFcstsT[j] = newsData.getOldForecast(targetSeries, endPerTarget.plus(j));
            revFcstsT[j] = newsData.getRevisedForecast(targetSeries, endPerTarget.plus(j));
            newFcstsT[j] = newsData.getNewForecast(targetSeries, endPerTarget.plus(j));
            
            if(!standardized){
                oldFcsts[j] = (oldFcstsT[j] * sdTarget) + meanTarget;
                revFcsts[j] = (revFcstsT[j] * sdTarget) + meanTarget;
                newFcsts[j] = (newFcstsT[j] * sdTarget) + meanTarget;
            }else{
                oldFcsts[j] = oldFcstsT[j];
                revFcsts[j] = revFcstsT[j];
                newFcsts[j] = newFcstsT[j];
            }  
        }           
        
        return DfmResultsNews.builder()
                .seriesIndex(index)
                .seriesPeriod(per)
                .seriesExpectedValueT(DoubleSeq.of(expT))
                .seriesExpectedValue(DoubleSeq.of(exp))
                .seriesObservedValueT(DoubleSeq.of(obsT))
                .seriesObservedValue(DoubleSeq.of(obs))
                .seriesNewsT(DoubleSeq.of(newsT))
                .seriesNews(DoubleSeq.of(news))
                .seriesWeightsT(sWeightsT)
                .seriesWeights(sWeights)
                .seriesImpactsT(sImpactsT)
                .seriesImpacts(sImpacts)
                .forecastsPeriods(fcstPeriods)
                .oldForecastsT(DoubleSeq.of(oldFcstsT))
                .oldForecasts(DoubleSeq.of(oldFcsts))
                .revisedForecastsT(DoubleSeq.of(revFcstsT))
                .revisedForecasts(DoubleSeq.of(revFcsts))
                .newForecastsT(DoubleSeq.of(newFcstsT))
                .newForecasts(DoubleSeq.of(newFcsts))
                .build();
    }
   

    
/*--------------------------------------------------------------------------
DEPRECATED METHODS (rjd3nowcasting v1): estimation of coefficients 
together with processing  
--------------------------------------------------------------------------*/
    
    public DfmResults estimate_PCA(DynamicFactorModel dfmModel, Matrix data, int freq, int[] start, boolean standardized, int nForecasts) {
        
        DfmResults.Builder builder = DfmResults.builder();
        
        TsInformationSet dfmData = prepareInput(data, freq, start, standardized, null, null, builder);
        
        PrincipalComponentsInitializer initializer = new PrincipalComponentsInitializer();
        DynamicFactorModel dfm = initializer.initialize(dfmModel, dfmData);

        DfmProcessor processor = DfmProcessor.builder().build();
        TsPeriod fLast = dfmData.getCurrentDomain().getEndPeriod().plus(nForecasts);
        TsInformationSet dfmDataExtended = dfmData.extendTo(fLast.start().toLocalDate());
        processor.process(dfm, dfmDataExtended);
        
        int nPeriod = dfmData.getCurrentDomain().getLength();
        int nPeriodExt = dfmDataExtended.getCurrentDomain().getLength();
        int nf = dfm.getNfactors();
        int neq =dfm.getMeasurementsCount();
        int blockLength = dfm.defaultSsfBlockLength();
        Matrix factors = getFactors(processor, nPeriodExt, nf, blockLength);
        Matrix factorsSdErr = getFactorsSdErr(processor, nPeriodExt, nf, blockLength);
        Matrix residuals = getResiduals(processor, dfmData, nPeriod, neq);
        Matrix residualsStandardized = getStandardizedResiduals(processor, dfmData, nPeriod, neq);
        // Not the same as GUI for the residuals... yet to investigate
        
        return builder
                .dfmData(dfmData)
                .dfm(dfm)
                .smoothedStates(processor.getSmoothingResults())
                .factors(factors)
                .factorsStdErr(factorsSdErr)
                .residuals(residuals)
                .residualsStandardized(residualsStandardized)
                .build();            
        // In case of PCA, nothing related to the likelihood
    }
    
    public DfmResults estimate_EM(DynamicFactorModel dfmModel, Matrix data, int freq, int[] start, boolean standardized, int nForecasts,
            boolean pcaInit, int maxIter, double eps) {
        
        DfmResults.Builder builder = DfmResults.builder();
        TsInformationSet dfmData = prepareInput(data, freq, start, standardized, null, null, builder);
        
        DynamicFactorModel dfm0;
        if (pcaInit) {
            PrincipalComponentsInitializer initializer = new PrincipalComponentsInitializer();
            DynamicFactorModel pcaModel = initializer.initialize(dfmModel, dfmData);
            if (pcaModel.isValid()) {
                dfm0 = pcaModel;
            }else{
                dfm0 = dfmModel;
            }
        }else{
            dfm0 = dfmModel;
        }
        
        DfmEM em = DfmEM.builder()
                .maxIter(maxIter)
                .precision(eps)
                .build();
        DynamicFactorModel dfm = em.initialize(dfm0, dfmData);
        
        DfmProcessor processor = DfmProcessor.builder().build();
        TsPeriod fLast = dfmData.getCurrentDomain().getEndPeriod().plus(nForecasts);
        TsInformationSet dfmDataExtended = dfmData.extendTo(fLast.start().toLocalDate());
        processor.process(dfm, dfmDataExtended);
        
        int nPeriod = dfmData.getCurrentDomain().getLength();
        int nPeriodExt = dfmDataExtended.getCurrentDomain().getLength();
        int nf = dfm.getNfactors();
        int neq =dfm.getMeasurementsCount();
        int blockLength = dfm.defaultSsfBlockLength();
        Matrix factors = getFactors(processor, nPeriodExt, nf, blockLength);
        Matrix factorsSdErr = getFactorsSdErr(processor, nPeriodExt, nf, blockLength);
        Matrix residuals = getResiduals(processor, dfmData, nPeriod, neq);
        Matrix residualsStandardized = getStandardizedResiduals(processor, dfmData, nPeriod, neq);
        
        return builder
                .dfmData(dfmData)
                .dfm(dfm)
                .smoothedStates(processor.getSmoothingResults())
                .factors(factors)
                .factorsStdErr(factorsSdErr)
                .residuals(residuals)
                .residualsStandardized(residualsStandardized)
                .logLikelihood(em.getFinalLogLikelihood())
                .build();
    }

    public DfmResults estimate_ML(DynamicFactorModel dfmModel, Matrix data, int freq, int[] start, boolean standardized, int nForecasts,
            boolean pcaInit, boolean emInit, int emInitmaxIter, double emInitEps, int maxIter, int maxBlockIter, int simplModelIter,
            boolean independentVARShocks, boolean mixedEstimation, double eps) {
            
        DfmResults.Builder builder = DfmResults.builder();
        TsInformationSet dfmData = prepareInput(data, freq, start, standardized, null, null, builder);
        TsPeriod fLast = dfmData.getCurrentDomain().getEndPeriod().plus(nForecasts);
        TsInformationSet dfmDataExtended = dfmData.extendTo(fLast.start().toLocalDate());
        
        DfmEstimator estimator = DfmEstimator.builder()
                .minimizer(ProxyMinimizer.builder(
                        LevenbergMarquardtMinimizer.builder()))
                .maxInitialIter(simplModelIter)
                .maxBlockIterations(maxBlockIter)
                .maxIterations(maxIter)
                .independentVarShocks(independentVARShocks)
                .mixed(mixedEstimation)
                .precision(eps)
                .build();
	
        DynamicFactorModel dfm0 = dfmModel;
        if (pcaInit) {
            PrincipalComponentsInitializer initializer = new PrincipalComponentsInitializer();
            dfm0 = initializer.initialize(dfmModel, dfmData);
        } 
        DfmKernel dfmK;
        if (emInit) {
            dfmK = DfmKernel.builder()
                    .initializer(DfmEM.builder()
                            .maxIter(emInitmaxIter)
                            .precision(emInitEps)
                            .build())
                    .estimator(estimator)
                    .processor(DfmProcessor.builder().build())
                    .build();
        } else {
            dfmK = DfmKernel.builder()
                    .estimator(estimator)
                    .processor(DfmProcessor.builder().build())
                    .build();
        }  
        dfmK.process(dfm0, dfmDataExtended);
        
        DynamicFactorModel dfm = dfmK.getEstimator().getEstimatedModel();
        DfmProcessor processor = dfmK.getProcessor();
        
        int nPeriod = dfmData.getCurrentDomain().getLength();
        int nPeriodExt = dfmDataExtended.getCurrentDomain().getLength();
        int nf = dfm.getNfactors();
        int neq = dfm.getMeasurementsCount();
        int blockLength = dfm.defaultSsfBlockLength();
        Matrix factors = getFactors(processor, nPeriodExt, nf, blockLength);
        Matrix factorsSdErr = getFactorsSdErr(processor, nPeriodExt, nf, blockLength);
        Matrix residuals = getResiduals(processor, dfmData, nPeriod, neq);
        Matrix residualsStandardized = getStandardizedResiduals(processor, dfmData, nPeriod, neq);

        return builder
                .dfmData(dfmData)
                .dfm(dfm)
                .smoothedStates(processor.getSmoothingResults())
                .factors(factors)
                .factorsStdErr(factorsSdErr)
                .residuals(residuals)
                .residualsStandardized(residualsStandardized)
                .gradient(dfmK.getEstimator().getGradient())
                .hessian(dfmK.getEstimator().getHessian())
                .likelihood(dfmK.getEstimator().getLikelihood())
                .logLikelihood(dfmK.getEstimator().getLikelihood().logLikelihood())
                .build();
    }
}


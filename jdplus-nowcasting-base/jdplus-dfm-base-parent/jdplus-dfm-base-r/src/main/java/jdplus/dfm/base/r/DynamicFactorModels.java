/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package jdplus.dfm.base.r;

import java.time.Month;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import jdplus.dfm.base.api.MeasurementType;
import jdplus.dfm.base.api.timeseries.TsInformationSet;
import jdplus.dfm.base.core.DfmEM;
import jdplus.dfm.base.core.DfmEstimates;
import jdplus.dfm.base.core.DfmEstimator;
import jdplus.dfm.base.core.DfmKernel;
import jdplus.dfm.base.core.DfmNews;
import jdplus.dfm.base.core.DfmProcessor;
import jdplus.dfm.base.core.DfmResults;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.dfm.base.core.DynamicFactorModel;
import jdplus.dfm.base.core.IDfmMeasurement;
import jdplus.dfm.base.core.MeasurementDescriptor;
import jdplus.dfm.base.core.PrincipalComponentsInitializer;
import jdplus.dfm.base.core.var.VarDescriptor;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.timeseries.TsCollection;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsDomain;
import jdplus.toolkit.base.api.timeseries.TsPeriod;
import jdplus.toolkit.base.core.data.DataBlock;
import jdplus.toolkit.base.core.data.DataBlockIterator;
import jdplus.toolkit.base.core.math.functions.levmar.LevenbergMarquardtMinimizer;
import jdplus.toolkit.base.core.math.functions.ssq.ProxyMinimizer;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.ssf.ISsfInitialization;
import jdplus.toolkit.base.core.stats.likelihood.Likelihood;
import jdplus.toolkit.base.r.timeseries.TsUtility;

/**
 *
 * @author LEMASSO
 */
@lombok.experimental.UtilityClass
public class DynamicFactorModels {

    private TsInformationSet prepareInput(Matrix data, int freq, int[] start, boolean standardized, DfmResults.Builder builder) {
        List<TsData> input = new ArrayList<>();
        double[] sampleMean = new double[data.getColumnsCount()];
        double[] sdDev = new double[data.getColumnsCount()];

        for (int j = 0; j < data.getColumnsCount(); ++j) {
            DoubleSeq sj = data.column(j);
            int nf = sj.count(Double::isFinite);
            double m = sj.averageWithMissing();
            double e2 = sj.ssqcWithMissing(m) / nf;
            double e = Math.sqrt(e2);

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
        System.out.println(dataT);
        
        if (builder != null) {
            builder.sampleMean(DoubleSeq.of(sampleMean))
                   .sampleStDev(DoubleSeq.of(sdDev))
                   .inputData(data)
                   .standardizedInput(standardized)
                   .transformedInputData(dataT);
        }        
        return tinput;
    }
    
    public Matrix getFactors(DfmProcessor processor, int nPeriodExtended, int nFactors, int mLags){
        FastMatrix factors = FastMatrix.make(nPeriodExtended, nFactors);       
        for (int j = 0; j < nFactors; ++j) {
            factors.column(j).add(processor.getSmoothingResults().getComponent((mLags+1)*j));  
        } 
        return factors;
    }
    
    public Matrix getFactorsSdErr(DfmProcessor processor, int nPeriodExtended, int nFactors, int mLags){
        FastMatrix factorsSdErr = FastMatrix.make(nPeriodExtended, nFactors);   
        for (int j = 0; j < nFactors; ++j) {
            factorsSdErr.column(j).add(processor.getSmoothingResults().getComponentVariance((mLags+1)*j).sqrt());
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
    
    // model from scratch
    public DynamicFactorModel model(int nfactors, int nlags, String[] factorType, Matrix factorLoaded, String varInit, double[] mVariance) {

        if (factorLoaded.getRowsCount() != factorType.length) {
            throw new IllegalArgumentException("The number of rows of factorLoaded should match the length of factor Type");
        }
        if (factorLoaded.getColumnsCount() != nfactors) {
            throw new IllegalArgumentException("The number of columns of factorLoaded should match nfactors");
        }

        int nc = nfactors * nlags;
        ISsfInitialization.Type varInitType = ISsfInitialization.Type.Unconditional;
        if (varInit.equals("Zero")) {
            varInitType = ISsfInitialization.Type.Zero;
        }
        VarDescriptor varDesc = new VarDescriptor(FastMatrix.make(nfactors, nc), FastMatrix.square(nfactors), varInitType);

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
                factorLoadedM.set(i, j, factorLoaded.get(i, j) == 0 ? Double.NaN : 1); // Because if not used, defined as NaN in class DynamicFactorModel
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
        if (!dfm.isValid()) {
            dfm = DynamicFactorModel.builder()
                    .measurements(mDescs)
                    .var(VarDescriptor.defaultVar(nfactors, nlags, varInitType))
                    .build();
        }
        return dfm;
    }
    
    // model with pre-initialized values of the parameters 
    public DynamicFactorModel model(int nfactors, int nlags, String[] factorType, Matrix factorLoaded, String varInit,
                                     Matrix varCoefficients, Matrix varInnovationsVariance, Matrix mCoefficients, double[] mIdiosyncraticErrVariance) {

        if (factorLoaded.getRowsCount() != factorType.length) {
            throw new IllegalArgumentException("The number of rows of factorLoaded should match the length of factor Type");
        }
        if (factorLoaded.getColumnsCount() != nfactors) {
            throw new IllegalArgumentException("The number of columns of factorLoaded should match nfactors");
        }
        
        int nc = nfactors * nlags;
        
        if(varCoefficients.getRowsCount() != nfactors || varCoefficients.getColumnsCount() != nc){
            throw new IllegalArgumentException("The dimension of the matrix varCoefficients is not consistent with the number of factors and/or lags.");
        }
        
        if(varInnovationsVariance.getRowsCount() != nfactors || !varInnovationsVariance.isSquare()){
            throw new IllegalArgumentException("The dimension of the matrix varInnovationsVariance is not consistent with the number of factors.");
        }
        
        if(mCoefficients.getRowsCount() != factorType.length || mCoefficients.getColumnsCount() !=  nfactors){
            throw new IllegalArgumentException("The dimension of the matrix parameters_factors is not consistent with the number of series (length of factorType) and/or the number of factors.");
        }
        
        if(mIdiosyncraticErrVariance.length != factorType.length){
            throw new IllegalArgumentException("The length of the vector parameters_factors_variance is not consistent with the number of series (length of factorType).");
        }
        
        ISsfInitialization.Type varInitType = ISsfInitialization.Type.Unconditional;
        if (varInit.equals("Zero")) {
            varInitType = ISsfInitialization.Type.Zero;
        }
        VarDescriptor varDesc = new VarDescriptor(varCoefficients, varInnovationsVariance, varInitType);

        FastMatrix mCoefficientsForm = FastMatrix.make(factorLoaded.getRowsCount(), factorLoaded.getColumnsCount());
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
                mCoefficientsForm.set(i, j, factorLoaded.get(i, j) == 0 ? Double.NaN : mCoefficients.get(i, j));
            }
            MeasurementDescriptor mdesc = MeasurementDescriptor.builder()
                    .type(factorTransf)
                    .coefficient(mCoefficientsForm.row(i))
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
                
        TsInformationSet dfmData = prepareInput(data, freq, start, standardized, null);
        
        PrincipalComponentsInitializer initializer = new PrincipalComponentsInitializer();
        DynamicFactorModel dfm = initializer.initialize(dfmModel, dfmData);
        
        DfmEstimates dfmE = DfmEstimates.builder()
                .dfm(dfm)
                .build();
        
        return dfmE; 
    }
      
    public DfmEstimates estimate_EM(DynamicFactorModel dfmModel, Matrix data, int freq, int[] start, boolean standardized, 
            boolean pcaInit, int maxIter, double eps) {
        
        TsInformationSet dfmData = prepareInput(data, freq, start, standardized, null);
        
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
            
        TsInformationSet dfmData = prepareInput(data, freq, start, standardized, null);
         
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
        
        TsInformationSet dfmData = prepareInput(data, freq, start, standardized, builder);
        
        DfmProcessor processor = DfmProcessor.builder().build();
        TsPeriod fLast = dfmData.getCurrentDomain().getEndPeriod().plus(nForecasts);
        TsInformationSet dfmDataExtended = dfmData.extendTo(fLast.start().toLocalDate());
        processor.process(dfmModel, dfmDataExtended);
        
        int nPeriod = dfmData.getCurrentDomain().getLength();
        int nPeriodExt = dfmDataExtended.getCurrentDomain().getLength();
        int nf = dfmModel.getNfactors();
        int neq =dfmModel.getMeasurementsCount();
        int mlags = dfmModel.measurementsLags();
        Matrix factors = getFactors(processor, nPeriodExt, nf, mlags);
        Matrix factorsSdErr = getFactorsSdErr(processor, nPeriodExt, nf, mlags);
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
        
    public boolean computeNews(DynamicFactorModel dfm, TsInformationSet oldData, TsInformationSet newData){
        
        DfmNews test = new DfmNews(dfm);
        test.process(oldData, newData);
        
        TsInformationSet dataOld = test.getOldInformationSet();
        TsInformationSet dataNew = test.getNewInformationSet();
        
        return true;
    }
    

    
    /*--------------------------------------------------------------------------
    DEPRECATED METHODS (rjd3nowcasting v1): estimation of coefficients 
    together with processing  
    --------------------------------------------------------------------------*/
    
    public DfmResults estimate_PCA(DynamicFactorModel dfmModel, Matrix data, int freq, int[] start, boolean standardized, int nForecasts) {
        
        DfmResults.Builder builder = DfmResults.builder();
        
        TsInformationSet dfmData = prepareInput(data, freq, start, standardized, builder);
        
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
        int mlags = dfm.measurementsLags();
        Matrix factors = getFactors(processor, nPeriodExt, nf, mlags);
        Matrix factorsSdErr = getFactorsSdErr(processor, nPeriodExt, nf, mlags);
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
        TsInformationSet dfmData = prepareInput(data, freq, start, standardized, builder);
        
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
        int mlags = dfm.measurementsLags();
        Matrix factors = getFactors(processor, nPeriodExt, nf, mlags);
        Matrix factorsSdErr = getFactorsSdErr(processor, nPeriodExt, nf, mlags);
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
        TsInformationSet dfmData = prepareInput(data, freq, start, standardized, builder);
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
        int mlags = dfm.measurementsLags();
        Matrix factors = getFactors(processor, nPeriodExt, nf, mlags);
        Matrix factorsSdErr = getFactorsSdErr(processor, nPeriodExt, nf, mlags);
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


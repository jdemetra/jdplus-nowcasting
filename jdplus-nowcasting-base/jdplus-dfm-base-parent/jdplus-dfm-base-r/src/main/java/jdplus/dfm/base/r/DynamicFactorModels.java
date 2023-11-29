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
import jdplus.dfm.base.api.NumericalProcessingSpec;
import jdplus.dfm.base.api.timeseries.TsInformationSet;
import jdplus.dfm.base.core.DfmEM;
import jdplus.dfm.base.core.DfmEstimator;
import jdplus.dfm.base.core.DfmKernel;
import jdplus.dfm.base.core.DfmProcessor;
import jdplus.dfm.base.core.DfmResults;
import jdplus.toolkit.base.api.information.GenericExplorable;
import jdplus.toolkit.base.api.information.InformationMapping;
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
        
        if (builder != null) {
            builder.sampleMean(DoubleSeq.of(sampleMean))
                   .sampleStDev(DoubleSeq.of(sdDev))
                   .inputData(data)
                   .standardizedInput(standardized)
                   .transformedInputData(tinput.generateMatrix(TsDomain.DEFAULT_EMPTY));
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
            boolean independantVAShocks, boolean mixedEstimation, double eps) {
            
        DfmResults.Builder builder = DfmResults.builder();
        TsInformationSet dfmData = prepareInput(data, freq, start, standardized, builder);
        TsPeriod fLast = dfmData.getCurrentDomain().getEndPeriod().plus(nForecasts);
        TsInformationSet dfmDataExtended = dfmData.extendTo(fLast.start().toLocalDate());
        
        DfmEstimator estimator = DfmEstimator.of(NumericalProcessingSpec.DEFAULT_ENABLED, dfmModel.getVar().getInitialization())
                .toBuilder()
                .maxIterations(maxIter)
                .maxBlockIterations(maxBlockIter)
                .maxInitialIter(simplModelIter)
                .independentVarShocks(independantVAShocks)
                .mixed(mixedEstimation)
                .precision(eps)
                .build();

        DfmKernel dfmK;
        if (pcaInit && emInit) {
            dfmK = DfmKernel.builder()
                    .initializer(DfmEM.builder()
                            .initializer(new PrincipalComponentsInitializer())
                            .maxIter(emInitmaxIter)
                            .precision(emInitEps)
                            .build())
                    .estimator(estimator)
                    .processor(DfmProcessor.builder().build())
                    .build();
        } else if (!pcaInit && emInit) {
            dfmK = DfmKernel.builder()
                    .initializer(DfmEM.builder()
                            .maxIter(emInitmaxIter)
                            .precision(emInitEps)
                            .build())
                    .estimator(estimator)
                    .processor(DfmProcessor.builder().build())
                    .build();
        } else if (pcaInit && !emInit) {
            dfmK = DfmKernel.builder()
                    .initializer(new PrincipalComponentsInitializer())
                    .estimator(estimator)
                    .processor(DfmProcessor.builder().build())
                    .build();
        } else {
            dfmK = DfmKernel.builder()
                    .estimator(estimator)
                    .processor(DfmProcessor.builder().build())
                    .build();
        }  
        dfmK.process(dfmModel, dfmDataExtended);
        
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
    
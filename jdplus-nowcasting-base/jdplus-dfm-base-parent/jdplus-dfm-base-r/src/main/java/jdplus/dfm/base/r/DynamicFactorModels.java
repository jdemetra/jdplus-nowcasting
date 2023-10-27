/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package jdplus.dfm.base.r;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import jdplus.dfm.base.api.MeasurementType;
import jdplus.dfm.base.api.NumericalProcessingSpec;
import jdplus.dfm.base.api.timeseries.TsInformationSet;
import jdplus.dfm.base.core.DfmEM;
import jdplus.dfm.base.core.DfmEstimator;
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
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.ssf.ISsfInitialization;
import jdplus.toolkit.base.r.timeseries.TsUtility;

/**
 *
 * @author LEMASSO
 */
@lombok.experimental.UtilityClass
public class DynamicFactorModels {

    public TsInformationSet prepareInput(Matrix data, int freq, int[] start, boolean standardized) {
        List<TsData> input = new ArrayList<>();
        for (int j = 0; j < data.getColumnsCount(); ++j) {
            DoubleSeq sj = data.column(j);
            if (!standardized) {
                int nf = sj.count(Double::isFinite);
                double m = sj.averageWithMissing();
                double e2 = sj.ssqcWithMissing(m) / nf;
                double e = Math.sqrt(e2);
                sj = sj.plus(-m);
                sj = sj.times(1 / e);
            }
            input.add(TsUtility.of(freq, start[0], start[1], sj.toArray()));
        }
        return new TsInformationSet(input);
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

    public DfmResults estimate_PCA(DynamicFactorModel dfmModel, Matrix data, int freq, int[] start, boolean standardized) {
        TsInformationSet dfmData = prepareInput(data, freq, start, standardized);
        PrincipalComponentsInitializer initializer = new PrincipalComponentsInitializer();
        DynamicFactorModel dfm = initializer.initialize(dfmModel, dfmData);

        DfmProcessor processor = DfmProcessor.builder().build();
        processor.process(dfm, dfmData);

        //System.out.println(dfm.getVar().getCoefficients());
        //System.out.println(dfm.getVar().getInnovationsVariance());
        return DfmResults.builder()
                .dfm(dfm)
                .smoothedStates(processor.getSmoothingResults())
                .build();
        // dans le cas du PCA, le DfmResults sera à moitié vide car rien par rapport au likelihood mais pas grave
    }

    public DfmResults estimate_EM(DynamicFactorModel dfmModel, Matrix data, int freq, int[] start, boolean standardized,
            boolean pcaInit, int maxIter, double eps) {

        TsInformationSet dfmData = prepareInput(data, freq, start, standardized);
        DynamicFactorModel dfm0 = dfmModel;
        if (pcaInit) {
            PrincipalComponentsInitializer initializer = new PrincipalComponentsInitializer();
            DynamicFactorModel pcaModel = initializer.initialize(dfmModel, dfmData);
            if (pcaModel.isValid()) {
                dfm0 = pcaModel;
            }
        }

        DfmEM em = DfmEM.builder()
                .maxIter(maxIter)
                .precision(eps)
                .build();
        DynamicFactorModel dfm = em.initialize(dfm0, dfmData);

        DfmProcessor processor = DfmProcessor.builder().build();
        processor.process(dfm, dfmData);

        System.out.println(dfm.getVar().getCoefficients());
        System.out.println(dfm.getVar().getInnovationsVariance());

        for (MeasurementDescriptor desc : dfm.getMeasurements()) {
            DoubleSeq coefficient = desc.getCoefficient();
            System.out.println(Arrays.toString(coefficient.toArray()));
            Double idiovar = desc.getVariance();
            System.out.println(idiovar);
        }

        return DfmResults.builder()
                .dfm(dfm)
                //.ll(em.getFinalLogLikelihood())
                .smoothedStates(processor.getSmoothingResults())
                .build();
        // je pourrai aussi retourner results de DynamicFactorModel grâce au delegate;

    }

    public DfmResults estimate_ML(DynamicFactorModel dfmModel, Matrix data, int freq, int[] start, boolean standardized,
            boolean pcaInit, boolean emInit, int emInitmaxIter, double emInitEps,
            int maxIter, int maxBlockIter, int simplModelIter, boolean independantVAShocks,
            boolean mixedEstimation, double eps) {

        TsInformationSet dfmData = prepareInput(data, freq, start, standardized);

        DynamicFactorModel dfm0 = dfmModel;
        if (pcaInit) {
            PrincipalComponentsInitializer initializer = new PrincipalComponentsInitializer();
            dfm0 = initializer.initialize(dfmModel, dfmData);
        }
        if (emInit) {
            DfmEM em = DfmEM.builder()
                    .maxIter(emInitmaxIter)
                    .precision(emInitEps)
                    .build();
            dfm0 = em.initialize(dfm0, dfmData);
        }

        DfmEstimator estimator = DfmEstimator.of(NumericalProcessingSpec.DEFAULT_ENABLED, dfm0.getVar().getInitialization())
                .toBuilder()
                .maxIterations(maxIter)
                .maxBlockIterations(maxBlockIter)
                .maxInitialIter(simplModelIter)
                .independentVarShocks(independantVAShocks)
                .mixed(mixedEstimation)
                .precision(eps)
                .build();
        estimator.estimate(dfm0, dfmData);
        DynamicFactorModel dfm = estimator.getEstimatedModel();

//        System.out.println(dfm.getVar().getCoefficients());
//        System.out.println(dfm.getVar().getInnovationsVariance());
//        
//        for (MeasurementDescriptor desc : dfm.getMeasurements()) {
//            DoubleSeq coefficient = desc.getCoefficient();
//            System.out.println(Arrays.toString(coefficient.toArray()));
//        }
        DfmProcessor processor = DfmProcessor.builder().build();
        processor.process(dfm, dfmData);

        return DfmResults.builder()
                .dfm(dfm)
                .likelihood(estimator.getLikelihood())
                .gradient(estimator.getGradient())
                .hessian(estimator.getHessian())
                .smoothedStates(processor.getSmoothingResults())
                .build();

// je pourrai aussi retourner results de DynamicFactorModel grâce au delegate;
    }

}

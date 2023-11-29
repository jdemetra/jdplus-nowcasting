
import jdplus.dfm.base.api.timeseries.TsInformationSet;
import jdplus.dfm.base.core.DfmEM;
import jdplus.dfm.base.core.DfmEstimator;
import jdplus.dfm.base.core.DfmKernel;
import jdplus.dfm.base.core.DfmResults;
import jdplus.dfm.base.core.DynamicFactorModel;
import jdplus.dfm.base.core.PrincipalComponentsInitializer;
import jdplus.dfm.base.r.DynamicFactorModels;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;

/**
 *
 * @author LEMASSO
 */
public class DynamicFactorModelsTest {

    public DynamicFactorModelsTest() {
    }

    public static void main(String[] args) {

        // Transformed data
        double[] PVI_FR = {-0.001, 0.003, 0, -0.003, 0.005, -0.001, -0.005, 0.005, 0.005, 0, -0.008, 0.001, 0.008, -0.003,
            -0.003, 0.004, -0.006, -0.001, -0.001, 0.008, -0.005, -0.001, 0.012, -0.007, 0.002, -0.004, 0.006, -0.001,
            0.009, -0.006, 0.003, 0.001, 0.003, 0.002, 0.001, 0.005, -0.016, 0.006, 0.004, -0.004, -0.005, 0.008, -0.001,
            0.001, -0.005, 0.005, -0.005, 0.007, 0.001, 0, -0.001, 0.002, 0.006, -0.011, 0.002, -0.005, 0.004, 0, -0.003,
            -0.008, 0.001, 0.009, -0.083, -0.101, 0.078, 0.051, 0.016, 0.002, 0.008, 0.011, -0.006, -0.002, 0.014, -0.02, 0.008,
            0, -0.001, 0, 0.002, 0, -0.005, 0.008, -0.007, 0, 0.009, -0.008, -0.001, -0.002, 0.001, 0.005, -0.007, 0.011, -0.004,
            -0.011, 0.008, 0.006, -0.006, 0.005, Double.NaN};
        double[] TURN_FR = {-0.003, 0.007, -0.002, 0.001, -0.004, 0.015, -0.013, -0.002, 0.005, -0.006, 0.005, -0.003, -0.001,
            -0.004, -0.001, 0.002, -0.001, 0.004, -0.003, 0.01, 0, -0.006, 0.015, 0.01, -0.014, 0.009, 0.018, -0.022, 0.012, -0.003,
            0.005, -0.007, 0.009, 0.004, 0.003, 0.005, -0.005, -0.009, 0.009, -0.001, -0.009, 0.014, -0.006, 0.009, -0.006, 0.019,
            -0.013, 0.006, -0.002, 0.008, -0.007, 0.005, 0.001, -0.016, 0.014, 0.003, -0.008, 0.007, -0.005, -0.004, -0.001, 0.004,
            -0.091, -0.109, 0.069, 0.07, 0.004, 0.01, 0.003, 0.011, 0, 0.004, 0.015, -0.018, 0.022, -0.003, -0.018, 0.033, -0.014,
            0.006, 0.007, 0.002, 0.012, 0.015, 0.008, 0, 0.01, -0.002, 0.005, 0.014, -0.01, 0.014, 0, -0.001, 0.007, 0.003, -0.002,
            Double.NaN, Double.NaN};
        double[] BS_FR = {-7.8, -7.3, -5.1, -5.9, -7, -7.1, -5.4, -5.7, -5.3, -3.9, -7.7, -4.8, -2.1, -0.3, -1.2, -3, -5.4, -6.7,
            -8.1, -8.1, -4.9, -4.6, -4.9, -3.3, -1.5, 0.4, -1.3, 1.1, 1.8, 1.2, 0.6, 1.3, 5.5, 3.8, 2.4, 4, 6.6, 4.7, 2.7, 5.2, 2.8,
            4.5, 2, 0.4, -0.4, -2.7, -1.3, -3.2, -2.7, -5, -3.9, -6.8, -1.8, -5.1, -7.9, -6.3, -7.2, -7, -7.5, -9.1, -2.8, -3, -8.6,
            -42.2, -26.9, -16.2, -13.4, -8.3, -9.9, -11.7, -16.1, -11.2, -10.1, -7.3, -7.1, -1.9, 2.1, 2.2, 6.8, 4.4, 2.3, 3.9, 5.5,
            4.9, 5.3, 7.5, 1.3, -0.2, 0.1, 0.2, -0.7, -4.6, -6.5, -6.8, -9.1, -8.5, -6.3, -6.3, -4.5};
        double[] PMI_FR = {0.4, 0, 1.2, -0.2, 0.2, 0.3, -0.1, -0.1, -0.3, 0.3, 0.5, 0.4, -0.9, -1.1, 0.4, 0.1, -0.2, 1.3, -0.8,
            -0.3, 0.9, 0.9, 0.2, 1.2, 0.3, 0.2, 0.8, 0.5, 0.3, 0.4, -0.8, 0.8, 0.7, 0.4, 1.6, 0.5, -1, -1, -2, -0.4, -0.7, -0.6, 0.2,
            -0.5, -1.4, -1.2, -0.2, -0.4, -0.9, -1.2, -1.8, 0.4, -0.2, -0.1, -1.1, 0.5, -1.3, 0.2, 1, -0.6, 1.6, 1.3, -4.7, -11.1, 6,
            8, 4.4, -0.1, 2, 1.1, -1, 1.4, -0.4, 3.1, 4.6, 0.4, 0.2, 0.3, -0.6, -1.4, -2.8, -0.3, 0.1, -0.4, 0.7, -0.5, -1.7, -1, -0.9,
            -2.5, -2.3, -0.1, -1.3, -2, 0.7, 0.7, 1, -0.3, -1.2};
        double[] GDP_FR = {Double.NaN, Double.NaN, 0.002, Double.NaN, Double.NaN, 0.001, Double.NaN, Double.NaN, 0.001, Double.NaN,
            Double.NaN, 0.001, Double.NaN, Double.NaN, 0.003, Double.NaN, Double.NaN, -0.002, Double.NaN, Double.NaN, 0.002, Double.NaN,
            Double.NaN, 0.003, Double.NaN, Double.NaN, 0.003, Double.NaN, Double.NaN, 0.004, Double.NaN, Double.NaN, 0.004, Double.NaN,
            Double.NaN, 0.003, Double.NaN, Double.NaN, 0, Double.NaN, Double.NaN, 0.002, Double.NaN, Double.NaN, 0.002, Double.NaN,
            Double.NaN, 0.002, Double.NaN, Double.NaN, 0.003, Double.NaN, Double.NaN, 0.003, Double.NaN, Double.NaN, 0, Double.NaN,
            Double.NaN, -0.002, Double.NaN, Double.NaN, -0.023, Double.NaN, Double.NaN, -0.061, Double.NaN, Double.NaN, 0.07, Double.NaN,
            Double.NaN, -0.003, Double.NaN, Double.NaN, 0, Double.NaN, Double.NaN, 0.004, Double.NaN, Double.NaN, 0.013, Double.NaN,
            Double.NaN, 0.002, Double.NaN, Double.NaN, 0, Double.NaN, Double.NaN, 0.002, Double.NaN, Double.NaN, 0.001, Double.NaN,
            Double.NaN, 0, Double.NaN, Double.NaN, 0};

        int[] start = {2015, 1};
        int nSeries = 5;
        int nRows = GDP_FR.length;
        double[][] seriesAll = new double[][]{PVI_FR, TURN_FR, BS_FR, PMI_FR, GDP_FR};
        FastMatrix data = FastMatrix.make(nRows, nSeries);
        int count = 0;
        while (count < seriesAll.length) {
            data.column(count).copyFrom(seriesAll[count], 0);
            count++;
        }

        // Model
        int nf = 2, nl = 2;
        String[] factorType = new String[]{"M", "M", "YoY", "M", "Q"};
        FastMatrix factorLoaded = FastMatrix.make(factorType.length, nf);
        double[] l1 = {1, 1};
        double[] l2 = {1, 0};
        double[] l3 = {0, 1};
        factorLoaded.row(0).copyFrom(l1, 0);
        factorLoaded.row(1).copyFrom(l3, 0);
        factorLoaded.row(2).copyFrom(l1, 0);
        factorLoaded.row(3).copyFrom(l2, 0);
        factorLoaded.row(4).copyFrom(l1, 0);
     
        DynamicFactorModel dfmInit = DynamicFactorModels.model(nf, nl, factorType, factorLoaded, "Unconditional", null);
        
        DfmResults dfm1 = DynamicFactorModels.estimate_PCA(dfmInit, data, 12, start, false, 12); // tested and same results as GUI
        //DfmResults dfm2 = DynamicFactorModels.estimate_EM(dfmInit, data, 12, start, false, 12, true, 100, 0.000000001); // tested and same results as GUI
        //DfmResults dfm3 = DynamicFactorModels.estimate_ML(dfmInit, data, 12, start, false, 12, false, false, 100, 0.000000001, 1000, 5, 15, false, true, 0.000000001);
        //Matrix test = dfm1.getHessian();
        //System.out.println(test.toString());
        

        DfmKernel kernel = DfmKernel.builder()
                .initializer(DfmEM.builder()
                        .initializer(new PrincipalComponentsInitializer())
                        .maxIter(100)
                        .build())
//                .estimator(DfmEstimator.builder()
//                        .build()
//                )
                .build();
//        kernel.process(dfmInit, DynamicFactorModels.prepareInput(data, 12, start, false));
//        long t0 = System.currentTimeMillis();
//        DfmResults dfm2 = DynamicFactorModels.estimate_EM(dfmInit, data, 12, start, false, 12, true, 100, 0.000000001); 
//        long t1 = System.currentTimeMillis();
//        System.out.println(t1 - t0);


        DfmResults dfm3 = DynamicFactorModels.estimate_ML(dfmInit, data, 12, start, false, 12, true, true,
                10, 0.0001, 500, 5, 10, false, true, 1e-9);
        System.out.println(dfm3.getDfm().getVar().getCoefficients());
        System.out.println(dfm3.getDfm().getVar().getInnovationsVariance());
        System.out.println(dfm3.getLikelihood().logLikelihood());
    }
}

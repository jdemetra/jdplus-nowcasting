
import jdplus.dfm.base.core.DfmEstimates;
import jdplus.dfm.base.core.DfmResults;
import java.util.ArrayList;
import java.util.List;
import jdplus.dfm.base.api.timeseries.TsInformationSet;
import jdplus.dfm.base.core.DfmNews;
import jdplus.dfm.base.core.DfmProcessor;
import jdplus.dfm.base.core.DfmResultsNews;
import jdplus.dfm.base.core.DynamicFactorModel;
import jdplus.dfm.base.r.DynamicFactorModels;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsPeriod;
import jdplus.toolkit.base.core.data.DataBlock;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;

/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
/**
 *
 * @author LEMASSO
 */
public class DynamicFactorModelsNewsTest {

    public DynamicFactorModelsNewsTest() {
    }

    public static void main(String[] args) {

        // Data vintages
        double[] PVI_FR_v1 = {-0.001, 0.003, 0, -0.003, 0.005, -0.001, -0.005, 0.005, 0.005, 0, -0.008, 0.001, 0.008, -0.003,
            -0.003, 0.004, -0.006, -0.001, -0.001, 0.008, -0.005, -0.001, 0.012, -0.007, 0.002, -0.004, 0.006, -0.001,
            0.009, -0.006, 0.003, 0.001, 0.003, 0.002, 0.001, 0.005, -0.016, 0.006, 0.004, -0.004, -0.005, 0.008, -0.001,
            0.001, -0.005, 0.005, -0.005, 0.007, 0.001, 0, -0.001, 0.002, 0.006, -0.011, 0.002, -0.005, 0.004, 0, -0.003,
            -0.008, 0.001, 0.009, -0.083, -0.101, 0.078, 0.051, 0.016, 0.002, 0.008, 0.011, -0.006, -0.002, 0.014, -0.02, 0.008,
            0, -0.001, 0, 0.002, 0, -0.005, 0.008, -0.007, 0, 0.009, -0.008, -0.001, -0.002, 0.001, 0.005, -0.007, 0.011, -0.004,
            -0.011, 0.008, 0.006, -0.006, Double.NaN};
        double[] PVI_FR_v2 = {-0.001, 0.003, 0, -0.003, 0.005, -0.001, -0.005, 0.005, 0.005, 0, -0.008, 0.001, 0.008, -0.003,
            -0.003, 0.004, -0.006, -0.001, -0.001, 0.008, -0.005, -0.001, 0.012, -0.007, 0.002, -0.004, 0.006, -0.001,
            0.009, -0.006, 0.003, 0.001, 0.003, 0.002, 0.001, 0.005, -0.016, 0.006, 0.004, -0.004, -0.005, 0.008, -0.001,
            0.001, -0.005, 0.005, -0.005, 0.007, 0.001, 0, -0.001, 0.002, 0.006, -0.011, 0.002, -0.005, 0.004, 0, -0.003,
            -0.008, 0.001, 0.009, -0.083, -0.101, 0.078, 0.051, 0.016, 0.002, 0.008, 0.011, -0.006, -0.002, 0.014, -0.02, 0.008,
            0, -0.001, 0, 0.002, 0, -0.005, 0.008, -0.007, 0, 0.009, -0.008, -0.001, -0.002, 0.001, 0.005, -0.007, 0.011, -0.004,
            -0.011, 0.008, 0.006, -0.006, 0.005, Double.NaN};
        double[] PVI_FR_v2bis = PVI_FR_v2;

        double[] TURN_FR_v1 = {-0.003, 0.007, -0.002, 0.001, -0.004, 0.015, -0.013, -0.002, 0.005, -0.006, 0.005, -0.003, -0.001,
            -0.004, -0.001, 0.002, -0.001, 0.004, -0.003, 0.01, 0, -0.006, 0.015, 0.01, -0.014, 0.009, 0.018, -0.022, 0.012, -0.003,
            0.005, -0.007, 0.009, 0.004, 0.003, 0.005, -0.005, -0.009, 0.009, -0.001, -0.009, 0.014, -0.006, 0.009, -0.006, 0.019,
            -0.013, 0.006, -0.002, 0.008, -0.007, 0.005, 0.001, -0.016, 0.014, 0.003, -0.008, 0.007, -0.005, -0.004, -0.001, 0.004,
            -0.091, -0.109, 0.069, 0.07, 0.004, 0.01, 0.003, 0.011, 0, 0.004, 0.015, -0.018, 0.022, -0.003, -0.018, 0.033, -0.014,
            0.006, 0.007, 0.002, 0.012, 0.015, 0.008, 0, 0.01, -0.002, 0.005, 0.014, -0.01, 0.014, 0, -0.001, 0.007, 0.003, -0.002,
            Double.NaN};
        double[] TURN_FR_v2 = {-0.003, 0.007, -0.002, 0.001, -0.004, 0.015, -0.013, -0.002, 0.005, -0.006, 0.005, -0.003, -0.001,
            -0.004, -0.001, 0.002, -0.001, 0.004, -0.003, 0.01, 0, -0.006, 0.015, 0.01, -0.014, 0.009, 0.018, -0.022, 0.012, -0.003,
            0.005, -0.007, 0.009, 0.004, 0.003, 0.005, -0.005, -0.009, 0.009, -0.001, -0.009, 0.014, -0.006, 0.009, -0.006, 0.019,
            -0.013, 0.006, -0.002, 0.008, -0.007, 0.005, 0.001, -0.016, 0.014, 0.003, -0.008, 0.007, -0.005, -0.004, -0.001, 0.004,
            -0.091, -0.109, 0.069, 0.07, 0.004, 0.01, 0.003, 0.011, 0, 0.004, 0.015, -0.018, 0.022, -0.003, -0.018, 0.033, -0.014,
            0.006, 0.007, 0.002, 0.012, 0.015, 0.008, 0, 0.01, -0.002, 0.005, 0.014, -0.01, 0.014, 0, -0.001, 0.007, 0.003, -0.002,
            Double.NaN, Double.NaN};
        double[] TURN_FR_v2bis = TURN_FR_v2;

        double[] BS_FR_v1 = {-7.8, -7.3, -5.1, -5.9, -7, -7.1, -5.4, -5.7, -5.3, -3.9, -7.7, -4.8, -2.1, -0.3, -1.2, -3, -5.4, -6.7,
            -8.1, -8.1, -4.9, -4.6, -4.9, -3.3, -1.5, 0.4, -1.3, 1.1, 1.8, 1.2, 0.6, 1.3, 5.5, 3.8, 2.4, 4, 6.6, 4.7, 2.7, 5.2, 2.8,
            4.5, 2, 0.4, -0.4, -2.7, -1.3, -3.2, -2.7, -5, -3.9, -6.8, -1.8, -5.1, -7.9, -6.3, -7.2, -7, -7.5, -9.1, -2.8, -3, -8.6,
            -42.2, -26.9, -16.2, -13.4, -8.3, -9.9, -11.7, -16.1, -11.2, -10.1, -7.3, -7.1, -1.9, 2.1, 2.2, 6.8, 4.4, 2.3, 3.9, 5.5,
            4.9, 5.3, 7.5, 1.3, -0.2, 0.1, 0.2, -0.7, -4.6, -6.5, -6.8, -9.1, -8.5, -6.3, -6.3};
        double[] BS_FR_v2 = {-7.8, -7.3, -5.1, -5.9, -7, -7.1, -5.4, -5.7, -5.3, -3.9, -7.7, -4.8, -2.1, -0.3, -1.2, -3, -5.4, -6.7,
            -8.1, -8.1, -4.9, -4.6, -4.9, -3.3, -1.5, 0.4, -1.3, 1.1, 1.8, 1.2, 0.6, 1.3, 5.5, 3.8, 2.4, 4, 6.6, 4.7, 2.7, 5.2, 2.8,
            4.5, 2, 0.4, -0.4, -2.7, -1.3, -3.2, -2.7, -5, -3.9, -6.8, -1.8, -5.1, -7.9, -6.3, -7.2, -7, -7.5, -9.1, -2.8, -3, -8.6,
            -42.2, -26.9, -16.2, -13.4, -8.3, -9.9, -11.7, -16.1, -11.2, -10.1, -7.3, -7.1, -1.9, 2.1, 2.2, 6.8, 4.4, 2.3, 3.9, 5.5,
            4.9, 5.3, 7.5, 1.3, -0.2, 0.1, 0.2, -0.7, -4.6, -6.5, -6.8, -9.1, -8.5, -6.3, -6.3, Double.NaN};
        double[] BS_FR_v2bis = {-7.8, -7.3, -5.1, -5.9, -7, -7.1, -5.4, -5.7, -5.3, -3.9, -7.7, -4.8, -2.1, -0.3, -1.2, -3, -5.4, -6.7,
            -8.1, -8.1, -4.9, -4.6, -4.9, -3.3, -1.5, 0.4, -1.3, 1.1, 1.8, 1.2, 0.6, 1.3, 5.5, 3.8, 2.4, 4, 6.6, 4.7, 2.7, 5.2, 2.8,
            4.5, 2, 0.4, -0.4, -2.7, -1.3, -3.2, -2.7, -5, -3.9, -6.8, -1.8, -5.1, -7.9, -6.3, -7.2, -7, -7.5, -9.1, -2.8, -3, -8.6,
            -42.2, -26.9, -16.2, -13.4, -8.3, -9.9, -11.7, -16.1, -11.2, -10.1, -7.3, -7.1, -1.9, 2.1, 2.2, 6.8, 4.4, 2.3, 3.9, 5.5,
            4.9, 5.3, 7.5, 1.3, -0.2, 0.1, 0.2, -0.7, -4.6, -6.5, -6.8, -9.1, -8.5, -6.3, -6.3, -4.5};

        double[] PMI_FR_v1 = {0.4, 0, 1.2, -0.2, 0.2, 0.3, -0.1, -0.1, -0.3, 0.3, 0.5, 0.4, -0.9, -1.1, 0.4, 0.1, -0.2, 1.3, -0.8,
            -0.3, 0.9, 0.9, 0.2, 1.2, 0.3, 0.2, 0.8, 0.5, 0.3, 0.4, -0.8, 0.8, 0.7, 0.4, 1.6, 0.5, -1, -1, -2, -0.4, -0.7, -0.6, 0.2,
            -0.5, -1.4, -1.2, -0.2, -0.4, -0.9, -1.2, -1.8, 0.4, -0.2, -0.1, -1.1, 0.5, -1.3, 0.2, 1, -0.6, 1.6, 1.3, -4.7, -11.1, 6,
            8, 4.4, -0.1, 2, 1.1, -1, 1.4, -0.4, 3.1, 4.6, 0.4, 0.2, 0.3, -0.6, -1.4, -2.8, -0.3, 0.1, -0.4, 0.7, -0.5, -1.7, -1, -0.9,
            -2.5, -2.3, -0.1, -1.3, -2, 0.7, 0.7, 1, -0.3};
        double[] PMI_FR_v2 = {0.4, 0, 1.2, -0.2, 0.2, 0.3, -0.1, -0.1, -0.3, 0.3, 0.5, 0.4, -0.9, -1.1, 0.4, 0.1, -0.2, 1.3, -0.8,
            -0.3, 0.9, 0.9, 0.2, 1.2, 0.3, 0.2, 0.8, 0.5, 0.3, 0.4, -0.8, 0.8, 0.7, 0.4, 1.6, 0.5, -1, -1, -2, -0.4, -0.7, -0.6, 0.2,
            -0.5, -1.4, -1.2, -0.2, -0.4, -0.9, -1.2, -1.8, 0.4, -0.2, -0.1, -1.1, 0.5, -1.3, 0.2, 1, -0.6, 1.6, 1.3, -4.7, -11.1, 6,
            8, 4.4, -0.1, 2, 1.1, -1, 1.4, -0.4, 3.1, 4.6, 0.4, 0.2, 0.3, -0.6, -1.4, -2.8, -0.3, 0.1, -0.4, 0.7, -0.5, -1.7, -1, -0.9,
            -2.5, -2.3, -0.1, -1.3, -2, 0.7, 0.7, 1, -0.3, Double.NaN};
        double[] PMI_FR_v2bis = {0.4, 0, 1.2, -0.2, 0.2, 0.3, -0.1, -0.1, -0.3, 0.3, 0.5, 0.4, -0.9, -1.1, 0.4, 0.1, -0.2, 1.3, -0.8,
            -0.3, 0.9, 0.9, 0.2, 1.2, 0.3, 0.2, 0.8, 0.5, 0.3, 0.4, -0.8, 0.8, 0.7, 0.4, 1.6, 0.5, -1, -1, -2, -0.4, -0.7, -0.6, 0.2,
            -0.5, -1.4, -1.2, -0.2, -0.4, -0.9, -1.2, -1.8, 0.4, -0.2, -0.1, -1.1, 0.5, -1.3, 0.2, 1, -0.6, 1.6, 1.3, -4.7, -11.1, 6,
            8, 4.4, -0.1, 2, 1.1, -1, 1.4, -0.4, 3.1, 4.6, 0.4, 0.2, 0.3, -0.6, -1.4, -2.8, -0.3, 0.1, -0.4, 0.7, -0.5, -1.7, -1, -0.9,
            -2.5, -2.3, -0.1, -1.3, -2, 0.7, 0.7, 1, -0.3, Double.NaN};

        double[] GDP_FR_v1 = {Double.NaN, Double.NaN, 0.002, Double.NaN, Double.NaN, 0.001, Double.NaN, Double.NaN, 0.001, Double.NaN,
            Double.NaN, 0.001, Double.NaN, Double.NaN, 0.003, Double.NaN, Double.NaN, -0.002, Double.NaN, Double.NaN, 0.002, Double.NaN,
            Double.NaN, 0.003, Double.NaN, Double.NaN, 0.003, Double.NaN, Double.NaN, 0.004, Double.NaN, Double.NaN, 0.004, Double.NaN,
            Double.NaN, 0.003, Double.NaN, Double.NaN, 0, Double.NaN, Double.NaN, 0.002, Double.NaN, Double.NaN, 0.002, Double.NaN,
            Double.NaN, 0.002, Double.NaN, Double.NaN, 0.003, Double.NaN, Double.NaN, 0.003, Double.NaN, Double.NaN, 0, Double.NaN,
            Double.NaN, -0.002, Double.NaN, Double.NaN, -0.023, Double.NaN, Double.NaN, -0.061, Double.NaN, Double.NaN, 0.07, Double.NaN,
            Double.NaN, -0.003, Double.NaN, Double.NaN, 0, Double.NaN, Double.NaN, 0.004, Double.NaN, Double.NaN, 0.013, Double.NaN,
            Double.NaN, 0.002, Double.NaN, Double.NaN, 0, Double.NaN, Double.NaN, 0.002, Double.NaN, Double.NaN, 0.001, Double.NaN,
            Double.NaN, 0, Double.NaN, Double.NaN};
        double[] GDP_FR_v2 = {Double.NaN, Double.NaN, 0.002, Double.NaN, Double.NaN, 0.001, Double.NaN, Double.NaN, 0.001, Double.NaN,
            Double.NaN, 0.001, Double.NaN, Double.NaN, 0.003, Double.NaN, Double.NaN, -0.002, Double.NaN, Double.NaN, 0.002, Double.NaN,
            Double.NaN, 0.003, Double.NaN, Double.NaN, 0.003, Double.NaN, Double.NaN, 0.004, Double.NaN, Double.NaN, 0.004, Double.NaN,
            Double.NaN, 0.003, Double.NaN, Double.NaN, 0, Double.NaN, Double.NaN, 0.002, Double.NaN, Double.NaN, 0.002, Double.NaN,
            Double.NaN, 0.002, Double.NaN, Double.NaN, 0.003, Double.NaN, Double.NaN, 0.003, Double.NaN, Double.NaN, 0, Double.NaN,
            Double.NaN, -0.002, Double.NaN, Double.NaN, -0.023, Double.NaN, Double.NaN, -0.061, Double.NaN, Double.NaN, 0.07, Double.NaN,
            Double.NaN, -0.003, Double.NaN, Double.NaN, 0, Double.NaN, Double.NaN, 0.004, Double.NaN, Double.NaN, 0.013, Double.NaN,
            Double.NaN, 0.002, Double.NaN, Double.NaN, 0, Double.NaN, Double.NaN, 0.002, Double.NaN, Double.NaN, 0.001, Double.NaN,
            Double.NaN, 0, Double.NaN, Double.NaN, Double.NaN};
        double[] GDP_FR_v2bis = GDP_FR_v2;

        int[] start = {2015, 1};
        int nSeries = 5;
        int nRows1 = GDP_FR_v1.length;
        int nRows2 = GDP_FR_v2.length;
        int nRows2bis = GDP_FR_v2bis.length;

        TsPeriod pstart = TsPeriod.monthly(2015, 1);
        double[][] seriesAllV1 = new double[][]{PVI_FR_v1, TURN_FR_v1, BS_FR_v1, PMI_FR_v1, GDP_FR_v1};
        double[][] seriesAllV2 = new double[][]{PVI_FR_v2, TURN_FR_v2, BS_FR_v2, PMI_FR_v2, GDP_FR_v2};
        double[][] seriesAllV2bis = new double[][]{PVI_FR_v2bis, TURN_FR_v2bis, BS_FR_v2bis, PMI_FR_v2bis, GDP_FR_v2bis};

        FastMatrix data1 = FastMatrix.make(nRows1, nSeries), data2 = FastMatrix.make(nRows2, nSeries), data2bis = FastMatrix.make(nRows2bis, nSeries);
        int count = 0;
        while (count < seriesAllV1.length) {
            data1.column(count).copyFrom(seriesAllV1[count], 0);
            data2.column(count).copyFrom(seriesAllV2[count], 0);
            data2bis.column(count).copyFrom(seriesAllV2bis[count], 0);
            count++;
        }

        // Model
        int nf = 2, nl = 2;
        String[] factorType = new String[]{"M", "M", "YoY", "M", "Q"};
        FastMatrix factorLoaded = FastMatrix.make(factorType.length, nf);
        double[] l1 = {1, 1};
        double[] l2 = {1, 0};
        double[] l3 = {1, 0};
        factorLoaded.row(0).copyFrom(l1, 0);
        factorLoaded.row(1).copyFrom(l3, 0);
        factorLoaded.row(2).copyFrom(l1, 0);
        factorLoaded.row(3).copyFrom(l2, 0);
        factorLoaded.row(4).copyFrom(l1, 0);

        DynamicFactorModel dfmInit = DynamicFactorModels.model(nf, nl, factorType, factorLoaded, "Unconditional", null);
        
//        DfmEstimates dfmEst = DynamicFactorModels.estimate_PCA(dfmInit, data1, 12, start, false);
        DfmEstimates dfmEst = DynamicFactorModels.estimate_EM(dfmInit, data1, 12, start, false, null, null, true, 100, 1e-09);
//        DfmEstimates dfmEst = DynamicFactorModels.estimate_ML(dfmInit, data1, 12, start, false, null, null, true, true, 100, 1e-09, 1000, 5, 15, false, true, 1e-09);
        
        DynamicFactorModel dfm = dfmEst.getDfm();    
        System.out.println(dfm);
        
        DfmResults rslt1 = DynamicFactorModels.process(dfm, data1, 12, start, false, null, null, 12);  
        DfmResults rslt2 = DynamicFactorModels.process(dfm, data2bis, 12, start, false, rslt1.getSampleMean().toArray(), rslt1.getSampleStDev().toArray(), 12); 
        
        
        // Tests differences in fcsts
        int nr = rslt2.getInputData().column(2).length();
        System.out.println("fcst1:");
        Matrix forecasts = rslt2.forecasts(1);
        System.out.println(forecasts.row(nr).get(2));
        
        DfmResultsNews test2 = DynamicFactorModels.computeNews(2, dfm, data1, data2bis, 12, start, false, 1);
        System.out.println("fcst2:");
        System.out.println(test2.getNewForecasts()); 
        
        
        // Tests class DfmNews
        List<TsData> ls1 = new ArrayList<>();
        List<TsData> ls2 = new ArrayList<>();
        List<TsData> ls2bis = new ArrayList<>();
        for (int i = 0; i < nSeries; ++i) {
            ls1.add(TsData.ofInternal(pstart, seriesAllV1[i]));
            ls2.add(TsData.ofInternal(pstart, seriesAllV2[i]));
            ls2bis.add(TsData.ofInternal(pstart, seriesAllV2bis[i]));
        }
        TsInformationSet ti1 = new TsInformationSet(ls1);
//        TsInformationSet ti2 = new TsInformationSet(ls2);
        TsInformationSet ti2bis = new TsInformationSet(ls2bis);
        
        System.out.println("");
        DfmNews test = new DfmNews(dfm);
        test.process(ti1, ti2bis);
        
        TsPeriod t = ls1.get(2).cleanExtremities().getEnd();
        System.out.println(test.revisions());
//        DoubleSeq rweights = test.weightsRevisions(4, t);
//        System.out.println(rweights);
        System.out.println(test.news());
        DoubleSeq weights = test.weights(2, t);
        System.out.println(weights);
        System.out.println(test.getOldForecast(4, t));
        System.out.println(test.getRevisedForecast(4, t));
        System.out.println(test.getNewForecast(4, t));
        System.out.println(test.revisions());
        
    }

}

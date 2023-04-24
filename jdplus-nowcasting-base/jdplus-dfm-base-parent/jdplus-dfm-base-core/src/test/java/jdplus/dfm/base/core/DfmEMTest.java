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

import java.io.File;
import java.io.IOException;

import tck.demetra.data.Data;
import tck.demetra.data.MatrixSerializer;
import java.util.ArrayList;
import java.util.List;
import jdplus.dfm.base.api.MeasurementType;
import jdplus.dfm.base.api.timeseries.TsInformationSet;
import jdplus.dfm.base.core.var.VarDescriptor;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsDomain;
import jdplus.toolkit.base.api.timeseries.TsPeriod;
import jdplus.toolkit.base.core.data.DataBlock;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.stats.DescriptiveStatistics;
import org.junit.jupiter.api.Test;

/**
 *
 * @author palatej
 */
public class DfmEMTest {

    static DynamicFactorModel dmodel;

    static final int N = 500;
    static final boolean stressTest = false;

    static TsInformationSet dfmdata;

    private static Matrix T, TVar, MVar, D, O, Z, M;
    private static FastMatrix dd, ddrnd;

//    // why static method? 
    private static void loadData() {   // DAVID: use a TsData[] to create a dfmdata object

        List<TsData> input = new ArrayList<>();
        TsPeriod start = TsPeriod.monthly(1980, 1);
        for (int i = 0; i < dd.getRowsCount(); ++i) {
            input.add(TsData.of(start, dd.row(i)));
        }
        dfmdata = new TsInformationSet(input); // input is TsData[]      
    }

    private static void loadDavidModel() throws IOException {
        File file = Data.copyToTempFile(DfmEMTest.class.getResource("/transition.csv"));
        T = MatrixSerializer.read(file, "\t|,");
//        System.out.println(T);
        file = Data.copyToTempFile(DfmEMTest.class.getResource("/tcovar.csv"));
        TVar = MatrixSerializer.read(file, "\t|,");
//        System.out.println(TVar);
        file = Data.copyToTempFile(DfmEMTest.class.getResource("/mcovar.csv"));
        MVar = MatrixSerializer.read(file, "\t|,");
//        System.out.println(MVar);
        file = Data.copyToTempFile(DfmEMTest.class.getResource("/data.csv"));
        D = MatrixSerializer.read(file, "\t|,");
        file = Data.copyToTempFile(DfmEMTest.class.getResource("/original.csv"));
        O = MatrixSerializer.read(file, "\t|,");
//         System.out.println(D);
        file = Data.copyToTempFile(DfmEMTest.class.getResource("/loadings.csv"));
        Z = MatrixSerializer.read(file, "\t|,");
        file = Data.copyToTempFile(DfmEMTest.class.getResource("/model.csv"));
        M = MatrixSerializer.read(file, "\t|,");
//        System.out.println(M);

        // build default the model
        //transition equation
        int nb = 3, nl = 4, c = 24, nc = nb * nl;
        FastMatrix C = FastMatrix.make(nb, nc);

        VarDescriptor var = new VarDescriptor(VarDescriptor.defaultCoefficients(nb, nl));
//
//        DynamicFactorModel.TransitionDescriptor tdesc = new DynamicFactorModel.TransitionDescriptor(nb, nl);
//        for (int i = 0; i < nb; ++i) {
//            DataBlock trow = TVar.row(i * 12).extract(0, nb, 12);
//            DataBlock prow = tdesc.covar.row(i);
//            prow.copy(trow);
//        }
//        dmodel.setTransition(tdesc);
//
        // measurement equation
        int nv = 0;
        for (int i = 0; i < M.getRowsCount(); ++i) {
            if (M.row(i).range(0, 3).anyMatch(x -> x != 0)) {
                ++nv;
            }
        }
        dd = FastMatrix.make(nv, O.getRowsCount() + 15);
        dd.set(Double.NaN);

        List<MeasurementDescriptor> mdescs = new ArrayList<>();

        for (int i = 0, j = 0; i < M.getRowsCount(); ++i) {
            if (M.row(i).range(0, 4).anyMatch(x -> x != 0)) {
                DataBlock row = dd.row(j);
                row.range(0, O.getRowsCount()).copy(O.column(j));
                DescriptiveStatistics stats = DescriptiveStatistics.of(O.column(j));
                double m = stats.getAverage();
                double e = stats.getStdev();
                row.sub(m);
                row.mul(1 / e);

                double[] q = new double[3];
                for (int k = 0; k < 3; ++k) {
                    if (M.get(i, k + 1) == 1) {
                        q[k] = Z.get(j, k * 12);
                    } else {
                        q[k] = Double.NaN;
                    }
                }
                MeasurementDescriptor desc = MeasurementDescriptor.builder()
                        .type(measurement((int) M.get(i, 0)))
                        .coefficient(DoubleSeq.of(q))
                        .variance(MVar.get(j, j))
                        .build();
                mdescs.add(desc);
                ++j;
            }
        }
        dmodel = new DynamicFactorModel(var, mdescs);
    }

    private static IDfmMeasurement measurement(int i) {
        switch (i) {
            case 1:
                return IDfmMeasurement.measurement(MeasurementType.M);
            case 2:
                return IDfmMeasurement.measurement(MeasurementType.Q);
            default:
                return IDfmMeasurement.measurement(MeasurementType.YoY);
        }
    }

    static {
        try//
        {
            loadDavidModel(); // first the model
            loadData();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    public DfmEMTest() {
    }

//    @Test
//   public void testinitCalc() {
//        DynamicFactorModel dmodelc = dmodel.clone();
//        dmodelc.normalize();dmodelc.setInitialization(VarSpec.Initialization.Zero);
//        PcInitializer initializer = new PcInitializer();
//       // initializer.initialize(dmodelc, dfmdata);
//
//
//      em = new DfmEM();
//      em.estimate(dmodelc,dfmdata);  // public, can be tested
//    //  em.estimate(dmodelc,dfmdata);  // public, can be tested
    @Test
    public void testInitialize() {
        PrincipalComponentsInitializer initializer = new PrincipalComponentsInitializer();
        TsDomain domain = dfmdata.getCurrentDomain().drop(120, 12);
        initializer.setEstimationDomain(domain);
        DynamicFactorModel model0 = initializer.initialize(dmodel, dfmdata);
        System.out.println(model0.getVar().getInnovationsVariance());
        System.out.println(model0.getVar().getCoefficients());
        for (MeasurementDescriptor desc : model0.getMeasurements()) {
            for (int i = 0; i < desc.getCoefficient().length(); ++i) {
                System.out.print(desc.getCoefficient(i));
                System.out.print('\t');
            }
            System.out.println(desc.getVariance());
        }
    }
}

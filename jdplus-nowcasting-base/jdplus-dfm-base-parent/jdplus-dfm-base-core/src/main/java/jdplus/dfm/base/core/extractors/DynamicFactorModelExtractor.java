/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package jdplus.dfm.base.core.extractors;

import jdplus.dfm.base.api.DfmDictionaries;
import jdplus.dfm.base.core.DynamicFactorModel;
import jdplus.toolkit.base.api.information.InformationExtractor;
import jdplus.toolkit.base.api.information.InformationMapping;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import nbbrd.service.ServiceProvider;

/**
 *
 * @author LEMASSO
 */
@ServiceProvider(InformationExtractor.class)
public class DynamicFactorModelExtractor extends InformationMapping<DynamicFactorModel> {
   
    public DynamicFactorModelExtractor() {
        set(DfmDictionaries.VAR_COEFFICIENTS, Matrix.class, source -> source.getVar().getCoefficients());
        set(DfmDictionaries.VAR_ERRORS_VARIANCE, Matrix.class, source -> source.getVar().getInnovationsVariance());
        set(DfmDictionaries.MEASUREMENT_COEFFICIENTS, Matrix.class, source -> {
            int nf = source.getNfactors();
            int nm = source.getMeasurementsCount();
            FastMatrix paramF = FastMatrix.make(nm, nf);
            for (int i = 0; i < nm; ++i) {
                paramF.row(i).copy(source.getMeasurements().get(i).getCoefficient());
            }
            return paramF;
        }
        );
        set(DfmDictionaries.MEASUREMENT_ERRORS_VARIANCE, double[].class, source -> {
            int nm = source.getMeasurementsCount();
            double[] varF = new double[nm];
            for (int i = 0; i < nm; ++i) {
                varF[i] = source.getMeasurements().get(i).getVariance();
            }
            return varF;
        }
        );
        
        set(DfmDictionaries.INITIALIZATION_TYPE, String.class, source -> source.getVar().getInitialization().toString());
        set(DfmDictionaries.FACTORS_TYPE, int[].class, source -> {
            int nm = source.getMeasurementsCount();
            int[] len = new int[nm];
            for (int i = 0; i < nm; ++i) {
                len[i] = source.getMeasurements().get(i).getType().getLength();
            }
            return len; 
        }
        );
    }
    
    @Override
    public Class<DynamicFactorModel> getSourceClass() {
        return DynamicFactorModel.class;
    }
}

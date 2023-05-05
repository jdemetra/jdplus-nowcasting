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

import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.core.data.DataBlock;
import static jdplus.toolkit.base.core.math.functions.IParametersDomain.PARAM;
import jdplus.toolkit.base.core.math.functions.IParametricMapping;
import jdplus.toolkit.base.core.math.functions.ParamValidation;



/**
 *
 * @author Jean
 */
public interface IDfmMapping extends IParametricMapping<DynamicFactorModel> {

    static final double EPS = 1e-6;

    DoubleSeq map(DynamicFactorModel m);
    
    @Override
    default double epsilon(DoubleSeq inparams, int idx) {
        double x = inparams.get(idx);
        return -x*EPS;
    }

    @Override
    default double lbound(int idx) {
        return -Double.MAX_VALUE;
    }

    @Override
    default double ubound(int idx) {
        return Double.MAX_VALUE;
    }

    @Override
    default ParamValidation validate(DataBlock ioparams) {
        return checkBoundaries(ioparams) ? ParamValidation.Valid : ParamValidation.Invalid;
    }

    @Override
    default String getDescription(int idx) {
        return PARAM + idx;
    }
}

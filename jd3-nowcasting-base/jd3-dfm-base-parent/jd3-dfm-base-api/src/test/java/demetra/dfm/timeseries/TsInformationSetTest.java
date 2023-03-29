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
package demetra.dfm.timeseries;

import demetra.data.AggregationType;
import demetra.math.matrices.Matrix;
import demetra.timeseries.TsData;
import demetra.timeseries.TsDomain;
import demetra.timeseries.TsPeriod;
import demetra.timeseries.TsUnit;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import static org.junit.Assert.*;
import tck.demetra.data.Data;

/**
 *
 * @author palatej
 */
public class TsInformationSetTest {
    
    public static final List<TsData> SERIES;
    
    static{
        SERIES=new ArrayList<>();
        SERIES.add(Data.TS_PROD);
        SERIES.add(Data.SP_IPI);
        SERIES.add(Data.TS_ABS_RETAIL);
        SERIES.add(Data.SP_IPI.aggregate(TsUnit.QUARTER, AggregationType.Average, true));
        SERIES.add(Data.TS_ABS_RETAIL2.aggregate(TsUnit.QUARTER, AggregationType.Sum, true));
        SERIES.add(Data.TS_PROD.aggregate(TsUnit.YEAR, AggregationType.Average, true).drop(0,3));
    }

    
    public TsInformationSetTest() {
    }

    @Test
    public void testMatrix() {
                
        TsInformationSet infoSet=new TsInformationSet(SERIES);
        Matrix M = infoSet.generateMatrix(TsDomain.of(TsPeriod.monthly(1975, 1), 300));
        M = infoSet.generateMatrix(TsDomain.of(TsPeriod.monthly(2000, 1), 1));
        M = infoSet.generateMatrix(TsDomain.of(TsPeriod.monthly(2000, 5), 100));
    }
    
}

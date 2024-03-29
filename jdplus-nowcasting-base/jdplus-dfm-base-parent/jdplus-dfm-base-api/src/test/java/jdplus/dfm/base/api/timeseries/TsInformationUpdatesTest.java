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
package jdplus.dfm.base.api.timeseries;

import java.util.ArrayList;
import java.util.List;
import jdplus.toolkit.base.api.data.AggregationType;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsDomain;
import jdplus.toolkit.base.api.timeseries.TsPeriod;
import jdplus.toolkit.base.api.timeseries.TsUnit;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;
import tck.demetra.data.Data;

/**
 *
 * @author palatej
 */
public class TsInformationUpdatesTest {

    public static final List<TsData> SERIES;
    public static final List<TsData> OSERIES;

    static {
        SERIES = new ArrayList<>();
        SERIES.add(Data.TS_PROD);
        SERIES.add(Data.SP_IPI);
        SERIES.add(Data.TS_ABS_RETAIL);
        SERIES.add(Data.SP_IPI.aggregate(TsUnit.QUARTER, AggregationType.Average, true));
        SERIES.add(Data.TS_ABS_RETAIL2.aggregate(TsUnit.QUARTER, AggregationType.Sum, true));
        SERIES.add(Data.TS_PROD.aggregate(TsUnit.YEAR, AggregationType.Average, true).drop(0, 3));

        OSERIES = new ArrayList<>();
        OSERIES.add(Data.TS_PROD.drop(0, 5));
        OSERIES.add(Data.SP_IPI.drop(0, 3));
        OSERIES.add(Data.TS_ABS_RETAIL.drop(0, 10));
        OSERIES.add(Data.SP_IPI.aggregate(TsUnit.QUARTER, AggregationType.Average, true).drop(0, 2));
        OSERIES.add(Data.TS_ABS_RETAIL2.aggregate(TsUnit.QUARTER, AggregationType.Sum, true).drop(0, 3));
        OSERIES.add(Data.TS_PROD.aggregate(TsUnit.YEAR, AggregationType.Average, true).drop(0, 4));
    }

    public TsInformationUpdatesTest() {
    }

    @org.junit.Test
    public void testCalendar() {

        TsInformationSet infoSet = new TsInformationSet(OSERIES);
        TsInformationSet ninfoSet = new TsInformationSet(SERIES);
        TsInformationUpdates updates = infoSet.updates(ninfoSet);
        assertTrue(updates.news().size() == 24);
        assertTrue(updates.revisions().isEmpty());

        TsPeriod first = updates.firstUpdate(12);
        TsPeriod last = updates.lastUpdate(12);
        assertTrue(last.isAfter(first));

        TsDomain udom = TsInformationUpdates.updatesDomain(1, updates.news());
        assertTrue(!udom.isEmpty());

        TsInformationSet rinfoSet = infoSet.revisedData(ninfoSet);
        for (int i = 0; i < rinfoSet.getSeriesCount(); ++i) {
            assertEquals(rinfoSet.series(i), infoSet.series(i));
        }
    }

}

/*******************************************************************************
*   Copyright 2014 Analog Devices, Inc.
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
********************************************************************************/

package com.analog.lyric.dimple.test.solvers.sumproduct;

import static java.util.Objects.*;
import static org.junit.Assert.*;

import org.junit.Test;

import com.analog.lyric.dimple.factorfunctions.And;
import com.analog.lyric.dimple.factorfunctions.Xor;
import com.analog.lyric.dimple.model.core.FactorGraph;
import com.analog.lyric.dimple.model.factors.Factor;
import com.analog.lyric.dimple.model.variables.Bit;
import com.analog.lyric.dimple.model.variables.Discrete;
import com.analog.lyric.dimple.options.BPOptions;
import com.analog.lyric.dimple.solvers.core.SolverBase;
import com.analog.lyric.dimple.solvers.gibbs.GibbsOptions;
import com.analog.lyric.dimple.solvers.minsum.MinSumSolver;
import com.analog.lyric.dimple.solvers.optimizedupdate.UpdateApproach;
import com.analog.lyric.dimple.solvers.sumproduct.Solver;
import com.analog.lyric.dimple.solvers.sumproduct.SumProductDiscrete;
import com.analog.lyric.dimple.solvers.sumproduct.SumProductSolver;
import com.analog.lyric.dimple.solvers.sumproduct.SumProductSolverGraph;
import com.analog.lyric.dimple.solvers.sumproduct.SumProductTableFactor;
import com.analog.lyric.dimple.solvers.sumproduct.sampledfactor.SampledFactor;
import com.analog.lyric.dimple.test.DimpleTestBase;

/**
 * 
 * @since 0.07
 * @author Christopher Barber
 */
public class TestSumProductOptions extends DimpleTestBase
{
	@SuppressWarnings({ "deprecation", "null" })
	@Test
	public void test()
	{
		// Test default values
		assertEquals(0.0, BPOptions.damping.defaultValue(), 0.0);
		assertTrue(BPOptions.nodeSpecificDamping.defaultValue().isEmpty());
		
		assertEquals(Integer.MAX_VALUE, (int)BPOptions.maxMessageSize.defaultValue());

		assertEquals(UpdateApproach.AUTOMATIC, BPOptions.updateApproach.defaultValue());
		assertEquals(1.0, BPOptions.automaticExecutionTimeScalingFactor.defaultValue(), 1.0e-9);
		assertEquals(10.0, BPOptions.automaticMemoryAllocationScalingFactor.defaultValue(), 1.0e-9);
		assertEquals(1.0, BPOptions.optimizedUpdateSparseThreshold.defaultValue(), 1.0e-9);
		
		final int nVars = 4;
		FactorGraph fg = new FactorGraph();
		Discrete[] vars = new Discrete[nVars];
		for (int i = 0; i < nVars; ++i)
		{
			vars[i] = new Bit();
		}
		Factor f1 = fg.addFactor(new Xor(), vars); // has custom factor
		Factor f2 = fg.addFactor(new And(), vars);
		
		// Check initial defaults
		SumProductSolverGraph sfg = requireNonNull(fg.setSolverFactory(new SumProductSolver()));
		assertEquals(0.0, sfg.getDamping(), 0.0);
		SumProductTableFactor sf1 = (SumProductTableFactor)requireNonNull(f1.getSolver());
		assertEquals(0, sf1.getK());
		assertEquals(0.0, sf1.getDamping(0), 0.0);
		SumProductTableFactor sf2 = (SumProductTableFactor)requireNonNull(f2.getSolver());
		assertEquals(0.0, sf2.getDamping(0), 0.0);
		assertEquals(0, sf2.getK());
		assertEquals(UpdateApproach.AUTOMATIC, sf1.getOptionOrDefault(BPOptions.updateApproach));
		
		assertEquals(SampledFactor.DEFAULT_BURN_IN_SCANS_PER_UPDATE, sfg.getSampledFactorBurnInScansPerUpdate());
		assertEquals(SampledFactor.DEFAULT_SAMPLES_PER_UPDATE, sfg.getSampledFactorSamplesPerUpdate());
		assertEquals(SampledFactor.DEFAULT_SCANS_PER_SAMPLE, sfg.getSampledFactorScansPerSample());
		
		assertNull(fg.setSolverFactory(null));
		
		// Set initial options on model
		fg.setOption(BPOptions.damping, .9);
		fg.setOption(BPOptions.maxMessageSize, 10);
		fg.setOption(GibbsOptions.burnInScans, 42); // will be overridden by default option in solver graph
		fg.setOption(GibbsOptions.scansPerSample, 23); // will be overridden by default option in solver graph
		fg.setOption(GibbsOptions.numSamples, 12); // will be overridden by default option in solver graph
		BPOptions.nodeSpecificDamping.set(f1, .4, .5, .6, .7);
		BPOptions.nodeSpecificDamping.set(f2, .3, .4, .5, .6);
		fg.setOption(BPOptions.updateApproach, UpdateApproach.AUTOMATIC);
		f2.setOption(BPOptions.updateApproach, UpdateApproach.NORMAL);
		
		// Test options that are updated on initialize()
		sfg = requireNonNull(fg.setSolverFactory(new SumProductSolver()));
		assertEquals(0.0, sfg.getDamping(), 0.0);
		assertEquals(0.0, sf1.getDamping(0), 0.0);
		sf1 = (SumProductTableFactor)requireNonNull(f1.getSolver());
		assertEquals(0, sf1.getK());
		sf2 = (SumProductTableFactor)requireNonNull(f2.getSolver());
		assertEquals(0, sf2.getK());
		assertEquals(SampledFactor.DEFAULT_BURN_IN_SCANS_PER_UPDATE, sfg.getSampledFactorBurnInScansPerUpdate());
		assertEquals(SampledFactor.DEFAULT_SAMPLES_PER_UPDATE, sfg.getSampledFactorSamplesPerUpdate());
		assertEquals(SampledFactor.DEFAULT_SCANS_PER_SAMPLE, sfg.getSampledFactorScansPerSample());
		assertEquals((Integer)SampledFactor.DEFAULT_SAMPLES_PER_UPDATE, sfg.getLocalOption(GibbsOptions.numSamples));
		assertEquals((Integer)SampledFactor.DEFAULT_SCANS_PER_SAMPLE, sfg.getLocalOption(GibbsOptions.scansPerSample));
		assertEquals((Integer)SampledFactor.DEFAULT_BURN_IN_SCANS_PER_UPDATE,
			sfg.getLocalOption(GibbsOptions.burnInScans));
		SumProductDiscrete sv1 = (SumProductDiscrete)vars[0].getSolver();
		assertEquals(0.0, sv1.getDamping(0), 0.0);
		
		sfg.initialize();
		assertEquals(.9, sfg.getDamping(), 0.0);
		assertEquals(.4, sf1.getDamping(0), 0.0);
		assertEquals(.5, sf1.getDamping(1), 0.0);
		assertEquals(.6, sf1.getDamping(2), 0.0);
		assertEquals(.7, sf1.getDamping(3), 0.0);
		assertEquals(10, sf1.getK());

		assertEquals(.3, sf2.getDamping(0), 0.0);
		assertEquals(.4, sf2.getDamping(1), 0.0);
		assertEquals(.5, sf2.getDamping(2), 0.0);
		assertEquals(.6, sf2.getDamping(3), 0.0);
		assertEquals(10, sf2.getK());
		
		assertEquals(.9, sv1.getDamping(0), 0.0);

		assertEquals(UpdateApproach.OPTIMIZED, sf1.getEffectiveUpdateApproach());
		assertEquals(UpdateApproach.NORMAL, sf2.getEffectiveUpdateApproach());
		
		// Test using set methods
		sfg.setDamping(.5);
		assertEquals(.5, sfg.getDamping(), 0.0);
		assertEquals(.5, requireNonNull(sfg.getLocalOption(BPOptions.damping)), 0.0);
		
		sf1.setK(3);
		assertEquals(3, sf1.getK());
		assertEquals((Integer)3, sf1.getLocalOption(BPOptions.maxMessageSize));
		
		sf1.setDamping(1, .23);
		assertEquals(.4, sf1.getDamping(0), 0.0);
		assertEquals(.23, sf1.getDamping(1), 0.0);
		assertEquals(.6, sf1.getDamping(2), 0.0);
		assertEquals(.7, sf1.getDamping(3), 0.0);
		assertArrayEquals(new double[] { .4,.23,.6,.7},
			BPOptions.nodeSpecificDamping.get(sf1).toPrimitiveArray(), 0.0);

		sfg.setSampledFactorSamplesPerUpdate(142);
		sfg.setSampledFactorScansPerSample(24);
		sfg.setSampledFactorBurnInScansPerUpdate(11);
		assertEquals(142, sfg.getSampledFactorSamplesPerUpdate());
		assertEquals(24, sfg.getSampledFactorScansPerSample());
		assertEquals(11, sfg.getSampledFactorBurnInScansPerUpdate());
		assertEquals((Integer)142, sfg.getLocalOption(GibbsOptions.numSamples));
		assertEquals((Integer)24, sfg.getLocalOption(GibbsOptions.scansPerSample));

		assertEquals(UpdateApproach.AUTOMATIC, sfg.getOption(BPOptions.updateApproach));

		assertNull(sf1.getLocalOption(BPOptions.updateApproach));
		sf1.setOption(BPOptions.updateApproach, UpdateApproach.OPTIMIZED);
		assertEquals(UpdateApproach.OPTIMIZED, sf1.getOption(BPOptions.updateApproach));
		assertEquals(UpdateApproach.OPTIMIZED, sf1.getLocalOption(BPOptions.updateApproach));
		sf1.setOption(BPOptions.updateApproach, UpdateApproach.NORMAL);
		assertEquals(UpdateApproach.NORMAL, sf1.getOption(BPOptions.updateApproach));
		assertEquals(UpdateApproach.NORMAL, sf1.getLocalOption(BPOptions.updateApproach));

		assertEquals((Integer)11, sfg.getLocalOption(GibbsOptions.burnInScans));
	}
	
	@Test
	public void testSolverEquality()
	{
		SolverBase<?> solver1 = new SumProductSolver();
		SolverBase<?> solver2 = new Solver();
		SolverBase<?> solver3 = new MinSumSolver();
		SolverBase<?> solver4 = new com.analog.lyric.dimple.solvers.gaussian.Solver();
		
		assertEquals(solver1, solver2);
		assertEquals(solver1.hashCode(), solver2.hashCode());
		assertNotEquals(solver1, solver3);
		assertNotEquals(solver1.hashCode(), solver3.hashCode());
		assertEquals(solver1, solver4);
		assertEquals(solver1.hashCode(), solver4.hashCode());
	}
}

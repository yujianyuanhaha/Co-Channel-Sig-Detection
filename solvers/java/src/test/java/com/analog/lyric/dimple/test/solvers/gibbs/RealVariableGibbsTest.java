/*******************************************************************************
*   Copyright 2012 Analog Devices, Inc.
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

package com.analog.lyric.dimple.test.solvers.gibbs;

import static java.util.Objects.*;
import static org.junit.Assert.*;

import org.junit.Test;

import com.analog.lyric.dimple.factorfunctions.MixedNormal;
import com.analog.lyric.dimple.factorfunctions.Normal;
import com.analog.lyric.dimple.model.core.FactorGraph;
import com.analog.lyric.dimple.model.variables.Discrete;
import com.analog.lyric.dimple.model.variables.Real;
import com.analog.lyric.dimple.solvers.gibbs.GibbsDiscrete;
import com.analog.lyric.dimple.solvers.gibbs.GibbsSolverGraph;
import com.analog.lyric.dimple.solvers.gibbs.GibbsReal;
import com.analog.lyric.dimple.test.DimpleTestBase;


public class RealVariableGibbsTest extends DimpleTestBase
{
	protected static boolean debugPrint = false;
	protected static boolean repeatable = true;
	
	@SuppressWarnings({ "deprecation", "null" })
	@Test
	public void basicTest1()
	{
		if (debugPrint) System.out.println("== basicTest1 ==");

		int numSamples = 10000;
		int updatesPerSample = 10;
		int burnInUpdates = 1000;

		FactorGraph graph = new FactorGraph();
		graph.setSolverFactory(new com.analog.lyric.dimple.solvers.gibbs.Solver());
		GibbsSolverGraph solver = (GibbsSolverGraph)graph.getSolver();
		solver.setNumSamples(numSamples);
		solver.setUpdatesPerSample(updatesPerSample);
		solver.setBurnInUpdates(burnInUpdates);
		
		double aPriorMean = 1;
		double aPriorSigma = 0.5;
		double aPriorR = 1/(aPriorSigma*aPriorSigma);
		double bPriorMean = -1;
		double bPriorSigma = 2.;
		double bPriorR = 1/(bPriorSigma*bPriorSigma);
		Real a = new Real();
		Real b = new Real();
		a.setInputObject(new Normal(aPriorMean,aPriorR));
		b.setInputObject(new Normal(bPriorMean, bPriorR));
		a.setName("a");
		b.setName("b");
		
		double abMean = 0;
		double abSigma = 1;
		double abR = 1/(abSigma*abSigma);
		graph.addFactor(new Normal(abMean,abR), a, b);
		
		GibbsReal sa = (GibbsReal)a.getSolver();
		GibbsReal sb = (GibbsReal)b.getSolver();
		
		sa.setProposalStandardDeviation(0.1);
		sb.setProposalStandardDeviation(0.1);

		
		if (repeatable) solver.setSeed(1);					// Make this repeatable
		solver.saveAllSamples();
		graph.solve();


		double[] aSamples = sa.getAllSamples();
		double[] bSamples = sb.getAllSamples();
		double aSum = 0;
		for (Object s : aSamples) aSum += (Double)s;
		double aMean = aSum/aSamples.length;
		if (debugPrint) System.out.println("aSampleMean: " + aMean);
		double bSum = 0;
		for (Object s : bSamples) bSum += (Double)s;
		double bMean = bSum/bSamples.length;
		if (debugPrint) System.out.println("bSampleMean: " + bMean);
		
		double aExpectedMean = (aPriorMean*aPriorR + abMean*abR)/(aPriorR + abR);
		double bExpectedMean = (bPriorMean*bPriorR + abMean*abR)/(bPriorR + abR);
		if (debugPrint) System.out.println("aExpectedMean: " + aExpectedMean);
		if (debugPrint) System.out.println("bExpectedMean: " + bExpectedMean);
		
		// Best should be the same as the mean in this case
		if (debugPrint) System.out.println("aBest: " + (Double)sa.getBestSample());
		if (debugPrint) System.out.println("bBest: " + (Double)sb.getBestSample());
		
		assertTrue(nearlyEquals(aMean,0.8050875226168582));
		assertTrue(nearlyEquals(bMean,-0.1921312702232493));
		assertTrue(nearlyEquals(sa.getBestSample(),0.8043550661413381));
		assertTrue(nearlyEquals(sb.getBestSample(),-0.20700427734616236));
	}
	
	
	@SuppressWarnings("deprecation")
	@Test
	public void basicTest2()
	{
		if (debugPrint) System.out.println("== basicTest2 ==");
		
		int numSamples = 10000;
		int updatesPerSample = 10;
		int burnInUpdates = 1000;

		FactorGraph graph = new FactorGraph();
		graph.setSolverFactory(new com.analog.lyric.dimple.solvers.gibbs.Solver());
		GibbsSolverGraph solver = requireNonNull((GibbsSolverGraph)graph.getSolver());
		solver.setNumSamples(numSamples);
		solver.setUpdatesPerSample(updatesPerSample);
		solver.setBurnInUpdates(burnInUpdates);
		
		double aPriorMean = 0;
		double aPriorSigma = 5;
		double aPriorR = 1/(aPriorSigma*aPriorSigma);
		double bProb1 = 0.6;
		double bProb0 = 1 - bProb1;
		Real a = new Real();
		a.setInputObject(new Normal(aPriorMean,aPriorR));
		Discrete b = new Discrete(0,1);
		b.setInput(bProb0, bProb1);
		a.setName("a");
		b.setName("b");
		
		double fMean0 = -1;
		double fSigma0 = 0.75;
		double fMean1 = 1;
		double fSigma1 = 0.75;
		double fR0 = 1/(fSigma0*fSigma0);
		double fR1 = 1/(fSigma1*fSigma1);
		graph.addFactor(new MixedNormal(fMean0, fR0, fMean1, fR1), a, b);
		
		GibbsReal sa = requireNonNull((GibbsReal)a.getSolver());
		GibbsDiscrete sb = requireNonNull((GibbsDiscrete)b.getSolver());
		
		sa.setProposalStandardDeviation(1.0);

		
		if (repeatable) solver.setSeed(1);					// Make this repeatable
		solver.saveAllSamples();
		graph.solve();


		double[] aSamples = sa.getAllSamples();
		Object[] bSamples = sb.getAllSamples();
		double aSum = 0;
		for (Object s : aSamples) aSum += (Double)s;
		double aMean = aSum/aSamples.length;
		if (debugPrint) System.out.println("aSampleMean: " + aMean);
		double bSum = 0;
		for (Object s : bSamples) bSum += (Integer)s;
		double bMean = bSum/bSamples.length;
		if (debugPrint) System.out.println("bSampleMean: " + bMean);
		
		if (debugPrint)
		{
			System.out.print("a = [");
			for (Object s : aSamples) System.out.print(s + " ");
			System.out.print("];\n");
		}

		
		double aExpectedMean = bProb0*(aPriorMean*aPriorR + fMean0*fR0)/(aPriorR + fR0) + bProb1*(aPriorMean*aPriorR + fMean1*fR1)/(aPriorR + fR1);
		if (debugPrint) System.out.println("aExpectedMean: " + aExpectedMean);
		if (debugPrint) System.out.println("bExpectedMean: " + bProb1);
		
		if (debugPrint) System.out.println("aBest: " + (Double)sa.getBestSample());
		if (debugPrint) System.out.println("bBest: " + sb.getBestSample());
		
		assertTrue(nearlyEquals(aMean,0.20867216566185906));
		assertTrue(nearlyEquals(bMean,0.6055));
		assertTrue(nearlyEquals(sa.getBestSample(),0.977986266650138));
		assertTrue((Integer)sb.getBestSample() == 1);
	}
	
	
	private static double TOLLERANCE = 1e-12;
	private boolean nearlyEquals(double a, double b)
	{
		double diff = a - b;
		if (diff > TOLLERANCE) return false;
		if (diff < -TOLLERANCE) return false;
		return true;
	}

}

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

package com.analog.lyric.dimple.solvers.core;

import static java.util.Objects.*;

import com.analog.lyric.collect.ArrayUtil;
import com.analog.lyric.dimple.model.factors.Factor;
import com.analog.lyric.dimple.solvers.interfaces.ISolverNode;
import com.analog.lyric.dimple.solvers.interfaces.ISolverVariable;

public abstract class STableFactorDoubleArray extends STableFactorBase
{
	protected double [][] _inputMsgs = ArrayUtil.EMPTY_DOUBLE_ARRAY_ARRAY;
	protected double [][] _outputMsgs = ArrayUtil.EMPTY_DOUBLE_ARRAY_ARRAY;

	public STableFactorDoubleArray(Factor factor)
	{
		super(factor);
	}

	@Override
	public void createMessages()
	{
		final Factor factor = _factor;
		int nVars = factor.getSiblingCount();
		
	    _inputMsgs = new double[nVars][];
	    _outputMsgs = new double[nVars][];
	    
	    for (int index = 0, end = nVars; index < end; index++)
	    {
	    	ISolverVariable is = requireNonNull(factor.getSibling(index).getSolver());
	    	Object [] messages = requireNonNull(is.createMessages(this));
	    	_outputMsgs[index] = (double[])messages[0];
	    	_inputMsgs[index] = (double[])messages[1];
	    }
	    
	}


	@Override
	public void resetEdgeMessages(int i)
	{
		
	}

	@Override
	public void moveMessages(ISolverNode other, int portNum, int otherPort)
	{

		STableFactorDoubleArray sother = (STableFactorDoubleArray)other;
	    _inputMsgs[portNum] = sother._inputMsgs[otherPort];
	    _outputMsgs[portNum] = sother._outputMsgs[otherPort];
	    
	}


	@Override
	public Object getInputMsg(int portIndex)
	{
		return _inputMsgs[portIndex];
	}

	@Override
	public Object getOutputMsg(int portIndex)
	{
		return _outputMsgs[portIndex];
	}

	@Override
	public void setInputMsgValues(int portIndex, Object obj)
	{
		double [] tmp = (double[])obj;
		for (int i = 0; i <tmp.length; i++)
			_inputMsgs[portIndex][i] = tmp[i];
	}
	
	@Override
	public void setOutputMsgValues(int portIndex, Object obj)
	{
		double [] tmp = (double[])obj;
		for (int i = 0; i <tmp.length; i++)
			_outputMsgs[portIndex][i] = tmp[i];
	}
	
}

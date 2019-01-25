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

package com.analog.lyric.dimple.model.variables;

import com.analog.lyric.util.misc.IVariableMapList;
import com.analog.lyric.util.misc.MapList;

public class VariableList extends MapList<Variable> implements IVariableMapList
{
	public VariableList()
	{
		super();
	}
	
	public VariableList(int initialCapacity)
	{
		super(initialCapacity);
	}

	/**
	 * Construct from an array.
	 * @since 0.07
	 */
	public VariableList(Variable [] vars)
	{
		super(vars.length);
		for (Variable v : vars)
			add(v);
	}

	/**
	 * Construct from an iterable.
	 * @since 0.07
	 */
	public VariableList(Iterable<Variable> vars)
	{
		super(vars);
	}

}

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

package com.analog.lyric.dimple.test.jsproxy;

import com.analog.lyric.dimple.jsproxy.IJSObject;
import com.analog.lyric.dimple.test.DimpleTestBase;

/**
 * 
 * @since 0.07
 * @author Christopher Barber
 */
public class JSTestBase extends DimpleTestBase
{
	final DimpleAppletTestState state = new DimpleAppletTestState();
	
	/**
	 * Creates a new JSObject for testing
	 * <p>
	 * @since 0.07
	 */
	IJSObject createJSObject(Object ... memberNamesAndValues)
	{
		FakeJSObject jsobj = new FakeJSObject();
		
		for (int i = 0; i < memberNamesAndValues.length; i+=2)
		{
			jsobj.setMember(memberNamesAndValues[i].toString(), memberNamesAndValues[i+1]);
		}
		
		return jsobj;
	}
}

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

import static com.analog.lyric.util.test.ExceptionTester.*;
import static org.junit.Assert.*;

import java.util.NoSuchElementException;

import org.junit.Test;

import com.analog.lyric.dimple.jsproxy.JSEnvironment;
import com.analog.lyric.dimple.jsproxy.JSOptionHolder;
import com.analog.lyric.dimple.options.BPOptions;
import com.analog.lyric.dimple.options.DimpleOptionRegistry;

/**
 * Test for {@link JSOptionHolder} class.
 * @since 0.07
 * @author Christopher Barber
 */
public class TestJSOptionHolder extends JSTestBase
{
	@Test
	public void test()
	{
		JSOptionHolder<?> holder = state.createGraph();
		
		String iterationsKey = BPOptions.iterations.qualifiedName();
		String dampingKey = BPOptions.damping.qualifiedName();
		
		expectThrow(NoSuchElementException.class, holder, "getOption", "does-not-exist");
		assertEquals(0.0, holder.getOption(dampingKey));
		assertEquals(1, holder.getOption(iterationsKey));
		assertFalse(holder.isOptionSet(dampingKey));
		
		holder.setOption(dampingKey, .9);
		assertTrue(holder.isOptionSet(dampingKey));
		assertEquals(.9, holder.getOption(dampingKey));
		assertEquals(.9, holder.getOption(BPOptions.damping));
		
		holder.setOption(iterationsKey, 42);
		assertEquals(42, holder.getOption(iterationsKey));
		assertTrue(holder.isOptionSet(iterationsKey));
		
		holder.unsetOption(iterationsKey);
		assertFalse(holder.isOptionSet(iterationsKey));
		assertEquals(1, holder.getOption(iterationsKey));
		assertEquals(.9, holder.getOption(dampingKey));
		
		holder.setOption(BPOptions.iterations, 23);
		assertEquals(23, holder.getOption(iterationsKey));
		
		holder.clearOptions();
		assertFalse(holder.isOptionSet(BPOptions.iterations));
		assertFalse(holder.isOptionSet(BPOptions.damping));
		
		JSEnvironment env = state.getEnvironment();
		assertEquals(state.applet, env.getApplet());
		env.setOption(iterationsKey, 123);
		assertEquals(123, holder.getOption(iterationsKey));

		DimpleOptionRegistry registry = env.getDelegate().optionRegistry();
		
		assertArrayEquals(registry.getAllMatching(".*").toArray(), (Object[])holder.getOptionKeysMatching(".*"));
		assertArrayEquals(registry.getAllMatching("damping").toArray(), (Object[])holder.getOptionKeysMatching("damping"));
	}
}

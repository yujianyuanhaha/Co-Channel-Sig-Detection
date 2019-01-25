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

package com.analog.lyric.dimple.jsproxy;

import com.analog.lyric.dimple.model.domains.RealDomain;

/**
 * Javascript API representation for a scalar, continuous real domain with optional bounds.
 * <p>
 * This wraps an underlying Dimple {@link RealDomain} object.
 * <p>
 * @since 0.07
 * @author Christopher Barber
 */
public class JSRealDomain extends JSDomain<RealDomain>
{
	JSRealDomain(JSDomainFactory factory, RealDomain domain)
	{
		super(factory, domain);
	}
	
	@Override
	public JSDomain.Type getDomainType()
	{
		return JSDomain.Type.REAL;
	}
	
	/**
	 * The lower bound of the domain.
	 * <p>
	 * All elements of the domain are greater than or equal to this value.
	 * @since 0.07
	 */
	public final double getLowerBound()
	{
		return getDelegate().getLowerBound();
	}

	/**
	 * The upper bound of the domain.
	 * <p>
	 * All elements of the domain are less than or equal to this value.
	 */
	public final double getUpperBound()
	{
		return getDelegate().getUpperBound();
	}
}

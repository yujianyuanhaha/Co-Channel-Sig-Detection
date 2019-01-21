function testOptions()
%testOptions Runs unit test of Dimple's option mechanism.
%
%   This only tests the overall option mechanism, not the behavior
%   of individual option settings.

env = DimpleEnvironment.active();
env.clearLocalOptions();

fg = FactorGraph();
v = Bit(2,2);
fg.addFactor('Xor',v);

% Make sure all options in dimpleOptions() are accepted.
optionNames = dimpleOptions();
for (i = 1:numel(optionNames))
    optionName = optionNames{i};
    val = fg.getOption(optionName);
    fg.setOption(optionName, val);
    fg.unsetOption(optionName);
end
assertTrue(isempty(getLocalOptions(fg)));

% Try a bogus key
try
    fg.getOption('Foo.bar');
    assertTrue(false, 'Expected an error');
catch err
end

assertEqual({1 1; 1 1}, v.getOption('BPOptions.iterations'));
v.setOption('BPOptions.iterations', 3);
assertEqual({3 3; 3 3}, v.getOption('BPOptions.iterations'));
v(1,1).setOption('BPOptions.iterations', 4);
v(2,:).setOption('BPOptions.iterations', 5);
assertEqual({4 3; 5 5}, v.getOption('BPOptions.iterations'));
assertEqual({3;5}, v(:,2).getOption('BPOptions.iterations'));
v(2,:).unsetOption('BPOptions.iterations');
assertEqual({4 3; 1 1}, v.getOption('BPOptions.iterations'));
v.clearLocalOptions();
assertEqual({1 1; 1 1}, v.getOption('BPOptions.iterations'));

env.setOption('BPOptions.iterations', 42);
assertEqual(42, env.getOption('BPOptions.iterations'));
assertEqual(42, v(1,1).getOption('BPOptions.iterations'));
v(2,2).setOption('BPOptions.iterations',12);
assertEqual({42 42; 42 12}, v.getOption('BPOptions.iterations'));

v.clearLocalOptions();
env.clearLocalOptions();
fg.clearLocalOptions();

assertTrue(isempty(fg.getLocalOptions()));
fg.setOption('BPOptions.iterations', 23);
fg.setOption('DimpleOptions.randomSeed', 42);
options = fg.getLocalOptions();
assertEqual({'BPOptions.iterations', 23; 'DimpleOptions.randomSeed', 42}, options);
fg.clearLocalOptions();
assertTrue(isempty(fg.getLocalOptions()));
fg.setOptions(options);
assertEqual({'BPOptions.iterations', 23; 'DimpleOptions.randomSeed', 42}, options);

v.setOptions('BPOptions.iterations', 21);
assertEqual(v(1,1).getOption('BPOptions.iterations'), 21);
v.setOptions('BPOptions.iterations', 27, 'DimpleOptions.randomSeed', 88);
assertEqual(v(1,2).getOption('BPOptions.iterations'), 27);
assertEqual(v(2,1).getOption('DimpleOptions.randomSeed'), 88);
v.setOptions({'BPOptions.iterations', 111, 'DimpleOptions.randomSeed', 34});
assertEqual(v(2,2).getOption('BPOptions.iterations'), 111);
assertEqual(v(1,1).getOption('DimpleOptions.randomSeed'), 34);

end


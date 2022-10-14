context( 'areBracketsBalanced')

expect_true(areBracketsBalanced("{}"), label="{}" )

expect_true(areBracketsBalanced("()"), label="()" )

expect_true(areBracketsBalanced("[]"), label="[]" )


expect_true(areBracketsBalanced("{{{
	[[[
	(((
	)))
	]]]
	}}}"), label="{{{[[[((()))]]]}}}" )


expect_true(!areBracketsBalanced("{{{[[((()))]]]}}}"), label="{{{[[((()))]]]}}}" )
expect_true(!areBracketsBalanced("{{[[[((()))]]]}}}"), label="{{[[[((()))]]]}}}" )
expect_true(!areBracketsBalanced("{{{[[[(()))]]]}}}"), label="{{{[[[(()))]]]}}}" )
expect_true(!areBracketsBalanced("{{{[[[((())]]]}}}"), label="{{{[[[((())]]]}}}" )
expect_true(!areBracketsBalanced("{{{[[[((()))]]}}}"), label="{{{[[[((()))]]}}}" )
expect_true(!areBracketsBalanced("{{{[[[((()))]]]}}"), label="{{{[[[((()))]]]}}" )
expect_true(!areBracketsBalanced("{{[[[((()))]]}}}"), label="{{[[[((()))]]}}}" )


expect_true(areBracketsBalanced("try this{ using that[(1,2),(1,2),(1,2)] and that[ (1,2) ] }"), label="try this1" )

expect_true(!areBracketsBalanced("try this{ using that[(1,2,(1,2),(1,2)] and that[ (1,2) ] }"), label="try this2" )

expect_true(!areBracketsBalanced("{(})"), label="{(})" )

expect_error( areBracketsBalanced( paste(rep("(", 255), rep(")",255), collapse="") ), "too many brackets -> max 500 !" )

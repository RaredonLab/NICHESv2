# rule set function test script

# setup
setwd("/Users/msbr/GitHub/NICHESMethods")
load_all()

# test rule set function application
test <- RuleSetFunction(query.set = NICHES::ncomms8866$Pair.Name,
                query.set.name = 'FANTOM5',
                special.separating.character = '_')

# view
View(test)

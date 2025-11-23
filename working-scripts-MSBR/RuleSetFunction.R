# Rule Set (Creates Feature Metdata)

RuleSetFunction <- function(query.set = NULL, # A vector of mechanisms, consisting of a sending part and receiving part, separated by a special character
                              query.set.name = 'FANTOM5',
                              special.separating.character = '_',
                            ligand.name = NULL,
                            receptor.name = NULL){

    # Initialize new metadata
    feature.metadata <- data.frame(QUERY.SET = query.set.name,
                                   MECHANISM = query.set,
                                   LIGAND = NA,
                                   RECEPTOR = NA,
                                   LIGAND.NAME = ligand.name,
                                   RECEPTOR.NAME = receptor.name,
                                   DEFINITION.LEVEL = NA,
                                  LIGAND.CLASS = NA,
                                  LIGAND.FAMILY = NA,
                                  RECEPTOR.CLASS = NA,
                                  RECEPTOR.FAMILY = NA)

    # Rule 1: MECHANISM is the feature label itself
    ## Performed above, upon feature.metadata initialization

    ## here is a change

    # Rule 2: LIGAND and RECEPTOR definition via special separating character
    feature.metadata$LIGAND <- data.frame(stringr::str_split_fixed(feature.metadata$MECHANISM,
                                               pattern = special.separating.character,
                                               n = 2))[,1]
    feature.metadata$RECEPTOR <- data.frame(stringr::str_split_fixed(feature.metadata$MECHANISM,
                                                                   pattern = special.separating.character,
                                                                   n = 2))[,2]

# Run functions
    feature.metadata <- DefineLigandFamily(feature.metadata)
    feature.metadata <- DefineReceptorFamily(feature.metadata)
    feature.metadata <- DefineLigandClass(feature.metadata)
    feature.metadata <- DefineReceptorClass(feature.metadata)


    ###### Format for output ####
    feature.metadata$LIGAND.CLASS <- factor(feature.metadata$LIGAND.CLASS)
    feature.metadata$RECEPTOR.CLASS <- factor(feature.metadata$RECEPTOR.CLASS)
    feature.metadata$LIGAND.FAMILY <- factor(feature.metadata$LIGAND.FAMILY)
    feature.metadata$RECEPTOR.FAMILY <- factor(feature.metadata$RECEPTOR.FAMILY)
    #### Return ####
    return(feature.metadata)
}

# test
ncomms8866.use <- NICHES::ncomms8866[NICHES::ncomms8866$Pair.Evidence %in% c('literature supported'),]
output <- RuleSetFunction(query.set = ncomms8866.use$Pair.Name,
                          ligand.name = ncomms8866.use$Ligand.Name,
                          receptor.name = ncomms8866.use$Receptor.Name)
to.do <- output[output$DEFINITION.LEVEL==0,]

View(to.do)
View(output)
View(ncomms8866.use)

# explore
table(output$LIGAND.CLASS)

temp <- output[output$LIGAND.CLASS == 'Growth Factor' & !is.na(output$LIGAND.CLASS),]
table(as.character(temp$LIGAND.FAMILY))

temp <- output[output$LIGAND.CLASS == 'Spatial Guidance' & !is.na(output$LIGAND.CLASS),]
table(as.character(temp$LIGAND.FAMILY))

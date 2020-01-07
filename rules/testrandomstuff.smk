######################################################################
######################################################################
#     Testing random stuff
######################################################################

rule testingrcode:
# extend common peaks between two replicates
	input:
		"testRdata.txt"
	output:
		"testRoutput.pdf"
	conda:
		"envs/renv.yml"
	params:
		color="red"
	script:
		"testRscript.R"

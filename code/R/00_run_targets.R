# Record the session for Docker
renv::init() 

# tar_manifest(fields = all_of("command"))
# tar_visnetwork()
targets::tar_make()
# inputs <- targets::tar_read(data)
# targets::tar_load(exp_tb)

## Take a snaphsot of all dependencies
renv::snapshot()


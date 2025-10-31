"""
    find_data_path(data_paths; local=false) -> String

Finds and returns the path to the data directory. If `local` is set to `true`, it checks for a directory named "data" in the current working directory and returns the path specified by `data_paths["local_data"]` if found.
If `local` is `false` or the local "data" directory is not found, it iterates through the `data_paths` dictionary and returns the first path that corresponds to an existing directory.

# Arguments
- `data_paths::Dict{String,String}`: A dictionary where keys are descriptive names of the paths and values are the actual paths to data directories.
- `local::Bool=false` (optional): A flag to indicate whether to look for a local "data" directory first.

# Returns
- `String`: The path to the data directory that was found.

# Throws
- `ArgumentError`: If `local` is `true` but no local "data" directory is found.

# Examples
```julia
data_paths = Dict("local_data" => "./data", "remote_data" => "/mnt/data")
path = find_data_path(data_paths, local=true)
```
This function is useful for applications that need to flexibly switch between local and remote data sources.
"""
function find_data_path(;local_data = false)
    data_storage_locations = Dict{String,String}(
    "external_data" => "C:\\Users\\ank10ki\\Desktop\\simulation_data_julia\\",
    "external_data_desk" => "C:\\Users\\Andreas\\Desktop\\simulation_data_julia\\",
    "knecht_data" => "C:\\Users\\Andreas\\Desktop\\simulation_data_julia_knecht\\")

    if local_data
        if isdir("data")
            return pwd()*"\\data\\"
        else
            throw(ArgumentError("No local data directory found"))
        end
    end
    for path in values(data_storage_locations)
        if isdir(path)
            return path 
        end
    end
end

"""
### `foldersize`

Calculates the size of a directory in megabytes. If `recursive` is set to `true`, it computes the total size by recursively summing the sizes of all files within the directory and its subdirectories. If `recursive` is `false`, it returns the size of the directory itself, which may not be meaningful as directories typically have a small, fixed size.

#### Parameters:
- `dirpath::String = pwd()`: The path to the directory whose size is to be calculated.Defaults to the current working directory.
- `recursive::Bool = true`: A flag indicating whether the size calculation should be recursive. Defaults to `true`.

#### Returns:
- `Float64`: The size of the directory in megabytes.

#### Example Usage:
```julia
size = foldersize("/path/to/directory", recursive=true)
```
This function is useful for monitoring disk usage by directories,
     especially in applications where managing data storage is critical.
"""
function foldersize(dirpath::String = pwd(); recursive::Bool = true)
	if !recursive
		return filesize(dirpath)
	else
		total = 0
		for (root,dirs,files) in walkdir(dirpath)
			total += filesize(root)
			for f in files
				file = joinpath(root,f)
				size = filesize(file)
				total += size
			end
		end
		return total/1024^3
	end
end



"""
### `length_vec`

Calculates and returns the Euclidean length (magnitude) of a 2-dimensional vector represented by a `Vector{Int64}`.

#### Parameters:
- `vec::Vector{Int64}`: A vector with two elements where `vec[1]` is the x-component and `vec[2]` is the y-component of the vector.

#### Returns:
- `Float64`: The Euclidean length of the vector.

#### Example Usage:
```julia
vector = [3, 4]
length = length_vec(vector) # Returns 5.0
``` 
"""
function length_vec(vec::Vector{Int64})
    len = sqrt(vec[1]^2+vec[2]^2)
    return len
end       
          

"""
@h methodname

Outputs documentations in jupyternotenbooks in VScode as markdown without bugs.

Example of how to use the `@h` macro:
```julia
@h res_scaling
```
Outputs documentations in jupyternotenbooks in VScode as markdown without bugs.
"""
macro h(x)
    quote
        display("text/markdown", @doc $(esc(x)))
    end    
end
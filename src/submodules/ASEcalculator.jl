module ASEcalculator

using ..ACEmd
using Sockets
using StaticArrays
# i-PI example https://github.com/i-pi/i-pi/blob/master/drivers/py/driver.py

export run_driver

const hdrlen = 12
const pos_type = typeof(SVector(1., 1., 1.)u"bohr") #should be bohr

function sendmsg(comm, message; nbytes=hdrlen)
    @info "Sending message" message
    if length(message) == nbytes
        final_message = message
    elseif length(message) < nbytes
        l = nbytes - length(message)
        final_message = message * repeat(' ', l)
    else
        error("Message is too long")
    end
    write(comm, message)
end


function recvmsg(comm, nbytes=hdrlen)
    raw_message = read(comm, nbytes)
    @assert length(raw_message) == nbytes  "recieved message was not correct"
    message = Char.(raw_message) |> String |> strip
    @info "Recieved message" message
    return message
end

function recvinit(comm)
    @info "Recieving INIT"
    bead_index = read(comm, Int32)
    nbytes = read(comm, Int32)
    raw_data = read(comm, nbytes)
    message = Char.(raw_data) |> String |> strip
    return (;
        :bead_index => bead_index,
        :message => message
    )
end


function recvposdata(comm)
    raw_cell = read(comm, 9*8)
    raw_icell = read(comm, 9*8) # drop this (inverce cell)
    natoms = read(comm, Int32)
    raw_pos = read(comm, 8*3*natoms)
    data_cell = reinterpret(pos_type, raw_cell)
    data_pos = reinterpret(pos_type, raw_pos)
    @info "Position data recieved"
    return (;
        :cell => Vector(data_cell),
        :positions => Vector(data_pos)  # clean type a little bit
    )
end


function sendforce(comm, e::Number, forces::AbstractArray, virial::AbstractMatrix)
    sendmsg(comm, "FORCEREADY")
    write(comm, ustrip(u"hartree", e) )
    write(comm, Int32( length(forces) ) )
    write(comm, usrip.(u"hartree/bohr", forces) )
    write(comm, usrip.(u"hartree", virial) )

    # Send single byte at end to make sure we are alive
    write(comm, one(Int32) )
    write(comm, zero(UInt8) )
end


function run_driver(address, pot::ACEmd.ACEpotential, init_structure; port=31415, unixpipe=false )
    if unixpipe
        comm = open("/tmp/ipi_"*address)
    else
        comm = connect(address, port)
    end
    has_init = false
    has_data = false
    data = nothing

    masses = atomic_mass(init_structure)
    symbols = atomic_symbol(init_structure)
    anumbers = atomic_number(init_structure)
    positions = position(init_structure)
    cell = bounding_box(init_structure)
    pbc =  boundary_conditions(init_structure)


    while true
        header = recvmsg(comm)

        if header == "STATUS"
            if !has_init
                sendmsg("NEEDINIT")
            elseif has_data
                sendmsg("HAVEDATA")
            else
                sendmsg("READY")
            end
        elseif header == "INIT"
            init = recvinit(comm)
            # we don't init anything for now
            has_init = true
        elseif header == "POSDATA"
            pos = recvposdata(comm)
            positions = pos[:positions]
            cell = pos[:cell]
            system = FastSystem(cell, pbc, positions, symbols, anumbers, masses)
            data = ace_energy_forces_virial(pot, system)
            has_data = true
        elseif header == "GETFORCE"
            sendforce(comm, data["energy"], data["forces"], data["virial"])
        elseif header == "EXIT"
            @info "Exiting calculator" 
            close(comm)
            break
        else
            close(comm)
            error("Message not recognised")
        end
        
    end
end

end
module IPIprotocol
# This is AtomsCalculators implementation for i-PI protocol.
# For reference see https://github.com/i-pi/i-pi/blob/master/drivers/py/driver.py

using AtomsBase
using AtomsCalculators
using Sockets
using StaticArrays
using Unitful
using UnitfulAtomic

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
    write(comm, final_message)
end


function recvmsg(comm, nbytes=hdrlen)
    raw_message = read(comm, nbytes)
    if length(raw_message) == 0
        @info "Server was probably closed and did not send EXIT"
        return "EXIT"
    end
    @assert length(raw_message) == nbytes "recieved message did not have correct lenght"
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


function sendforce(comm, e::Number, forces::AbstractVector, virial::AbstractMatrix)
    etype = (eltype ∘ eltype)(forces)
    f_tmp = reinterpret(reshape, etype, forces)
    sendmsg(comm, "FORCEREADY")
    write(comm, ustrip(u"hartree", e) )
    write(comm, Int32( length(forces) ) )
    write(comm, ustrip.(u"hartree/bohr", f_tmp) )
    write(comm, ustrip.(u"hartree", virial) )

    # Send single byte at end to make sure we are alive
    write(comm, one(Int32) )
    write(comm, zero(UInt8) )
end

"""
    run_driver(address, potential, init_structure; port=31415, unixsocket=false )

Connect I-PI driver to server at given `address`. Use kword `port` (default 31415) to
specify port. If kword `unixsocket` is true, `address` is understood to be the name of the socket
and `port` option is ignored.
"""
function run_driver(address, potential, init_structure; port=31415, unixsocket=false )
    if unixsocket
        comm = connect("/tmp/ipi_"*address)
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
                sendmsg(comm, "NEEDINIT")
            elseif has_data
                sendmsg(comm, "HAVEDATA")
            else
                sendmsg(comm, "READY")
            end
        elseif header == "INIT"
            init = recvinit(comm)
            # we don't init anything for now
            has_init = true
            has_data = false
        elseif header == "POSDATA"
            pos = recvposdata(comm)
            positions = pos[:positions]
            cell = pos[:cell]
            @assert length(symbols) == length(positions) "received amount of position data does no match the atomic symbol data"
            system = FastSystem(cell, pbc, positions, symbols, anumbers, masses)
            data = AtomsCalculators.energy_forces_virial(system, potential)
            has_data = true
        elseif header == "GETFORCE"
            sendforce(comm, data[:energy], data[:forces], data[:virial])
            has_data = false
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
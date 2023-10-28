using Sockets
# i-PI example https://github.com/i-pi/i-pi/blob/master/drivers/py/driver.py


const hdrlen = 12
const pos_type = typeof(SVector(1., 1., 1.)u"pm") #should be bohr

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
    raw_icell = read(comm, 9*8) # drop this
    natoms = read(comm, Int32)
    raw_pos = read(comm, 8*3*natoms)
    data_cell = reinterpret(pos_type, raw_cell)
    data_pos = reinterpret(pos_type, raw_pos)
    @info "Position data recieved"
    return (;
        :cell => Vector(data_cell)
        :positions => Vector(data_pos)  # clean type a little bit
    )
end


function sendforce(comm, e::Number, forces::AbstractArray, virial::AbstractMatrix)
    sendmsg(comm, "FORCEREADY")
    write(comm, uconvert(u"eV", e) )
    write(comm, Int32( length(forces) ) )
    write(comm, uconvert(u"eV/pm", forces) )
    write(comm, uconvert(u"eV*pm", virial) )

    # Send single byte at end to make sure we are alive
    write(comm, one(Inte32) )
    write(comm, zero(UInt8) )
end


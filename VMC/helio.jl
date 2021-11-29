# Importamos las paqueterías a utilizar
using Random
using Plots
using LinearAlgebra

#Parámetros para la simulación
number_walkers = 50
steps = 60000
therm_step = 600
therm = 10000
trials1 = 24
dR = 1.0

#Parámetro a variar para la función de onda de prueba
A1 = 1.25
dA1 = 0.025

# Almacenamiento de puntos para gráfica
X = Vector(0:(trials1-1))
X = A1 .+ X .*dA1
Y = Vector{Float64}()


#Funciones de generación de puntos aleatorios
DIST = 2
function ran_pos()
    return [DIST .* (2 .* rand(3) .-1), DIST .* (2 .* rand(3) .-1)]
end

function aleatorio(n)
    rand(n).*2.0 .-1.0
end

function distancia(R::Vector{Vector{Float64}})
    suma = [0.0,0.0,0.0]
    suma = R[2] .- R[1]
    norm(suma)
end

function dot(R::Vector{Vector{Float64}}, i) 
    s = [0.0,0.0,0.0]
    s = R[i].*(R[1] .- R[2])
    suma = s[1] + s[2] + s[3]
    suma/(norm(R[i])*distancia(R))
end

#Función para generar nuevo paso
function paso(R::Vector{Vector{Float64}})
    return [R[1] + aleatorio(3).*dR, R[2] + aleatorio(3).*dR]
end

# Distribución de probabilidad
function prob(R::Vector{Vector{Float64}}, A::Vector{Float64})
    r1 = norm(R[1])
    r2 = norm(R[2])
    return (exp.(- A *(r1 + r2))).^2 
end

#Energía local
function energia(R::Vector{Vector{Float64}}, A::Float64)
    r = distancia(R)
    r1 = norm(R[1])
    r2 = norm(R[2])
    V = (2/r1 + 2/r2 - 1/r)
    -1 * A^2 + A / r1 + A /r2 - V
end

#Obtener el valor mínimo
minimo = [0.0,0.0]

function metropolis(R::Vector{Vector{Float64}}, A::Float64)
    points = []
    for i in 0:(steps-1)
        R_nuevo = paso(R)
        p = 0 
        try
            p = prob(R_nuevo, A)/prob(R,A)
        catch DivideError
            
        end
        if p >= rand()
            R = R_nuevo
        end
        if i>therm
            if (i-therm) % therm_step == 0
                push!(points, R)
            end
        end
    end
    return points
end

function int(points::Vector{Any}, A::Float64)
    suma = energia.(points, A)./length(points)
end 


function helio()
    A = Vector(0:(trials1-1))
    A = A1 .+ A .* dA1

    R = Vector{Vector{Vector{Float64}}}(undef, (trials1*number_walkers))
    R .= ran_pos.()

    for i in 1:(number_walkers-1)
        a = Vector(0:(trials1-1))
        a = A1 .+ a .* dA1
        A = [A;a]
    end
    
    m = metropolis.(R,A)
    s = int.(m,A)
    e = sum.(s)
    for i in 1:length(A)
        if e[i] < minimo[1]
            minimo[1] = e[i]
            minimo[2] = A[i]
        end
    end
    println(minimo)
end

helio()
#using Plots
using LinearAlgebra
screen = Matrix{Int64}(undef,640,480)
res = (480,640)

X = Array{Float64}([1,0,0])
Y = Array{Float64}([0,1,0])
Z = Array{Float64}([0,0,1])

struct ray 
  origin::Array{Float64}
  #p2::Array{Float64}
  dir::Array{Float64}
  i::Int
  j::Int
end
#ray(o,p) = ray(o,p,(p.-o)/norm(p.-o);i=0,j=0)
struct plane
  point::Array{Float64}
  normal::Array{Float64}
end
struct sphere
  origin::Array{Float64}
  radius::Float64
  source::Bool
end

struct camera
  pos::Array
  focus::Float64
  dirz::Array
  dirx::Array
  diry::Array
  rays::Array{ray}
  screen::Matrix
  shscreen::Matrix
end

function camera(pos,focus,pixd,dirz=Z,dirx=X,diry=Y)
  screen = 1 ./ zeros((res[1],res[2]))
  shscreen = 1 ./ zeros((res[1],res[2]))
  rays = Array{ray}([])
  scr_mid = pos + focus*dirz
  for i in 1:res[1]
    for j in 1:res[2]
      p = pos + focus*dirz + (j-res[2]/2)*pixd*dirx+ (i-res[1]/2)*pixd*diry
      #println(p)
      push!(rays,ray(pos,(p.-pos) ./ norm(p.-pos),i,j))
    end
  end
  camera(pos,focus,dirz,dirx,diry,rays,screen,shscreen)
end

point_on_ray(d::Float64,r::ray) = (r.origin.+d*r.dir)

#objects = [plane([0,-1,0],[0,1,0]), sphere([-.5,-0.5,0],3,false), sphere([10,10,0],1,true)] 
objects = [plane([0,-1,0],[0,1,0]), sphere([-.5,-0.5,0],3,false), sphere([0,0,-2],1,true)] 

function solve(r::ray,s::sphere,dmax::Float64=30,shadow::Bool=false)
  o = r.origin
  u = r.dir
  c = s.origin
  rad = s.radius
  nabla = dot(u,(o.-c))^2 - (norm(o.-c)^2-rad^2) 
  d = 0.0
  sh = 0.0
  if nabla < 0.0
     d = dmax
  else
    d1 = -dot(u,(o.-c))+sqrt(nabla )
    d2 = -dot(u,(o.-c))-sqrt(nabla)
    if d1<0 && d2<0
      d = dmax
    elseif d1>0 && d2<0
      d = d1
    elseif d1<0 && d2>0
      d = d2
    elseif d1>0 && d2>0
      d = min(d1,d2)
    end
  end
  po = point_on_ray(d,r)
  normal = (po-s.origin)/norm(po-s.origin)
  newdir = r.dir - 2*dot(r.dir,normal)*normal 
  newray = ray(po,objects[end].origin,0,0) #shadow ray
  dist_to_source = sqrt(sum((objects[end].origin .- po).^2))
  if shadow
    #for o in [obj for obj in objects if !(obj in [s,objects[end]])]
    for o in [objects[1], ]
      if 0<solve(newray,o,dmax,false)[1]<dist_to_source
        sh = 1.0
      end
    end
  end
  d<0 ? (return dmax,sh) : (return clamp(d,0,dmax),sh)
end

function solve(r::ray,p::plane,dmax::Float64=30,shadow::Bool=false)
  l=r.dir
  l0 = r.origin
  n = p.normal
  p0 = p.point
  d = dot((p0.-l0),n)/dot(l,n)
  po = point_on_ray(d,r)
  sh = 0
  #newdir = r.dir - 2*dot(r.dir,p.normal)*p.normal 
  #newray = ray(po,newdir,0,0)
  newray = ray(po,objects[3].origin,0,0) #shadow ray
  dist_to_source = sqrt(sum((objects[3].origin .- po).^2))
  if shadow
    #for o in [obj for obj in objects if !(obj in [p,objects[end]])]
    for o in [objects[2],]
      if 0<solve(newray,o,dmax,false)[1]<dist_to_source
        sh = 1.0
      end
    end
  end
  d<0 ? (return dmax,sh) : (return clamp(d,0,dmax),sh)
end

function draw_on_screen(cam,objects)
  for obj in objects[1:end]
    for r in cam.rays
      d,sh = solve(r,obj,30.0,true)
      if d < (cam.screen[r.i , r.j] )  
	 cam.screen[r.i , r.j ] = d 
	 cam.shscreen[r.i,r.j ] = sh
      end
    end
  end
end

#@time ray_intersect(cam.rays[45],objects[1],objects)
objects = [plane([0,-1,0],[0,1,0]), sphere([0,-0.0,0],1,false), sphere([7,0,0],0.5,true)] 
cam = camera(-3Z,0.1, 5e-4)
#cam.screen .=  1 ./ zeros((res[1],res[2]))
draw_on_screen(cam,objects)
heatmap(1 ./(cam.screen).^2,color=:greys)
heatmap(cam.shscreen)
heatmap(1 ./(cam.screen).^2 .- 0.1cam.shscreen,color=:greys)
#heatmap(256 .- cam.screen,color=:greys)

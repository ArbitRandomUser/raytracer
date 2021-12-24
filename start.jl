using Plots
using FileIO,ImageMagick,Colors,FixedPointNumbers
using LinearAlgebra
res = (480,640)
X = Array{Float64}([1,0,0])
Y = Array{Float64}([0,1,0])
Z = Array{Float64}([0,0,1])
O = Array{Float64}([0,0,0])

abstract type object end
struct ray 
  origin::Array{Float64}
  dir::Array{Float64}
  i::Int
  j::Int
end
struct plane <:object
  point::Array{Float64}
  normal::Array{Float64}
  shadow::Bool
  refl::Bool
end
struct sphere <:object
  origin::Array{Float64}
  radius::Float64
  shadow::Bool
  refl::Bool
end

struct camera
  pos::Array
  focus::Float64
  dirz::Array
  dirx::Array
  diry::Array
  rays::Array{ray}
  dscreen::Matrix
  shscreen::Matrix
  d2sscreen::Matrix
end
function camera(pos,focus,pixd,dirz=Z,dirx=X,diry=Y)
  screen = 1 ./ zeros((res[1],res[2]))
  shscreen = 1 ./ zeros((res[1],res[2]))
  d2sscreen = 1 ./ zeros((res[1],res[2]))
  rays = Array{ray}([])
  scr_mid = pos + focus*dirz
  for i in 1:res[1]
    for j in 1:res[2]
      p = pos + focus*dirz + (j-res[2]/2)*pixd*dirx+ (i-res[1]/2)*pixd*diry
      #println(p)
      push!(rays,ray(pos,(p.-pos) ./ norm(p.-pos),i,j))
    end
  end
  camera(pos,focus,dirz,dirx,diry,rays,screen,shscreen,d2sscreen)
end
function reset_camera(cam::camera)
  cam.screen .= 1 ./ zeros((res[1],res[2]))
  cam.shscreen .= 1 ./ zeros((res[1],res[2]))
end

function normalat(p,o::sphere)
  (p.-o.origin)/norm(p.-o.origin)
end
function normalat(p,o::plane)
  o.normal
end

function reflection(r,p,o,dmax=30.0)
  newdir = r.dir -2*dot(r.dir,normalat(p,o))*normalat(p,o) 
  newray = ray(p,newdir,0,0)
  ref_d = dmax
  for ob in objects[1:end-1]
    if ob!=o
      for d in solve(newray,ob,dmax)
	if 0<d<=dmax && d<ref_d
	  ref_d = d
	end
      end
    end
  end
  return ref_d
end

function shadow(r,po,o,dmax=30.0)
   newdir = (objects[end] .-po);newdir = newdir/norm(newdir)
   newray = ray(po,newdir,0,0) #shadow ray
   d2s = sqrt(sum((objects[end] .- po).^2))
   for ob in objects[1:end-1]# if !(obj in [p,objects[end]])]
     d2 = solve(newray,ob,dmax)
     if any(1e-9 .< d2 .< d2s)
       return 1.0 , d2s
     end
   end
   return 0.0 , d2s
end

point_on_ray(d::Float64,r::ray) = (r.origin.+d*r.dir)
function first_intersect(r::ray,o::object,dmax::Float64=30.0,)
  darr =  solve(r,o,dmax)
  ret_d = dmax 
  sh = 0.0 
  d2s = dmax
  for d in darr
    if 0<d<=dmax && d<ret_d 
      ret_d = d
    end
  end
  po = point_on_ray(ret_d,r)
  #if o.shadow #&& ret_d != dmax
  sh,d2s = shadow(r,po,o,dmax)
  if o.refl
    ret_d = ret_d + 0.3*reflection(r,po,o,dmax)
  end
  return (ret_d,sh,d2s)
end


function solve(r::ray,s::sphere,dmax::Float64=30)#,shadow::Bool=false)
  o,u,c,rad = r.origin,r.dir,s.origin,s.radius
  nabla = dot(u,(o.-c))^2 - (norm(o.-c)^2-rad^2) 
  if nabla <0
    return [dmax,dmax]
  else
    d1 = -dot(u,(o.-c))+sqrt(nabla )
    d2 = -dot(u,(o.-c))-sqrt(nabla)
    return [d1,d2]
  end
end

function solve(r::ray,p::plane,dmax::Float64=30)#shadow::Bool=false)
  l,l0,n,p0 = r.dir,r.origin,p.normal,p.point
  d = dot((p0.-l0),n)/dot(l,n)
  return [d,]
end

function draw_on_screen(cam,objects)
  for obj in objects[1:end-1]
     Threads.@threads for r in cam.rays
      d,sh,d2s = first_intersect(r,obj,30.0)
      if d < (cam.dscreen[r.i , r.j] )  
	 cam.dscreen[r.i , r.j ] = d 
	 cam.shscreen[r.i,r.j ] = sh !=1.0 ? sh : sh 
	 cam.d2sscreen[r.i,r.j] = d2s
      end
    end
  end
end
invsqr(arr) =  1 ./ arr.^2

lightsource = [2,2.0,-0.0]
floorwall = plane(-1.1Y,+Y,true,false)
leftwall = plane(-1.1X,+X,true,false)
rightwall = plane(+X,-X,true,false)
ball1 = sphere(O,1.0,true,true)
ball2 = sphere(O+0.56X-0.81Z-0.73Y,0.27,true,false)
objects = [floorwall,leftwall,ball1,ball2,lightsource]
cam = camera(-2.5Z,0.1, 5e-4)
draw_on_screen(cam,objects)
#imgarr = (1 ./(1*cam.dscreen+0.3cam.shscreen).^2 )[end:-1:begin,:]
imgarr = ( invsqr(cam.d2sscreen+0.5cam.shscreen) )
imgarr = (imgarr)[end:-1:begin,:]
imgarr = imgarr/(invsqr(1.5))
save("test.png",Array{RGB{N0f8}}( imgarr))


#let
#nseconds = 5 
#nframes = 24*nseconds
#gifarr = zeros((res[1],res[2],nframes))
#mov = 10.0/nframes 
#count = 0
#prev_t = time()
#for i in 1:nframes
#  lightsource .+= [0.0,0.0,mov]
#  reset_camera(cam)
#  draw_on_screen(cam,objects)
#  imgarr = (1 ./(cam.screen+0.3cam.shscreen).^2 )[end:-1:begin,:]
#  gifarr[:,:,i] .= imgarr 
#  println(i," frames rendered")
#  println(time() - prev_t)
#  prev_t = time()
#end
#save("test.gif",Array{RGB{N0f8}}(gifarr),fps=24)
#end

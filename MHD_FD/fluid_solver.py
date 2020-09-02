import numpy as np
from solver import solver_2D


class fluid_solver_2D(solver_2D):

    def linear_advect_explicit_2d(self, f):
        return self.dt*self.bd_1st_2d(f, self.C, self.C)

    def linear_advect_implicit_2d(self, f, XX, YY):
        x = XX - (self.dt/self.dx * self.C)
        y = YY - (self.dt/self.dy * self.C)

        x = np.where(x < 0.5, 0.5, x)
        y = np.where(y < 0.5, 0.5, y)

        x = np.where(x > (self.Nx-2)+0.5, (self.Nx-2)+0.5, x)
        y = np.where(x > (self.Nx-2)+0.5, (self.Ny-2)+0.5, y)

        io = x.astype(int)
        jo = y.astype(int)
        i1 = io + 1
        j1 = jo + 1

        s1 = x - io
        t1 = y - jo
        so = 1 - s1
        to = 1 - t1

        return (so*(to*f[io, jo]+t1*f[io, j1])+s1*(to*f[i1, jo]+t1*f[i1, j1]))

    def linear_advect_implicit_periodic_2d(self, f, XX, YY):

        x = XX - (self.dt/self.dy * C)
        y = YY - (self.dt/self.dy * C)

        x = x % (self.Nx - 2)
        y = y % (self.Nx - 2)

        io = x.astype(int)
        jo = y.astype(int)

        i1 = io + 1
        j1 = jo + 1

        s1 = x - io
        t1 = y - jo

        s0 = 1 - s1
        to = 1 - t1

        return (so*(to*f[io, jo]+t1*f[io, j1])+s1*(to*f[i1, jo]+t1*f[i1, j1]))

    def nonlinear_advect_explicit_2d(self, f, fx, fy):
        return self.dt*self.bd_1st_2d(f, fx[1:-1, 1:-1], fy[1:-1, 1:-1])

    def nonlinear_advect_implicit_2d(self, f, fx, fy, XX, YY):
        x = XX - (self.dt/self.dx * fx[1:-1, 1:-1])
        y = YY - (self.dt/self.dx * fy[1:-1, 1:-1])

        x = np.where(x < 0.5, 0.5, x)
        y = np.where(y < 0.5, 0.5, y)

        x = np.where(x > (self.Nx-2) + 0.5, (self.Nx-2)+0.5, x)
        y = np.where(x > (self.Ny-2) + 0.5, (self.Ny-2)+0.5, y)

        io = x.astype(int)
        jo = y.astype(int)

        i1 - io + 1
        j1 = jo + 1

        s1 = x - io
        t1 = y - jo
        so = 1 - s1
        to = 1 - t1

        return (so*(to*f[io, jo]+t1*f[io, j1])+s1*(t0*f[i1, jo]+t1*f[i1, j1]))

    def nonlinear_advect_implicit_periodic_2d(self, f, fx, fy, XX, YY):

        x = XX - (self.dt/self.dx * fx[1:-1, 1:-1])
        y = YY - (self.dt/self.dy * fy[1:-1, 1:-1])

        x = x % (self.Nx - 2)
        y = y % (self.Ny - 2)

        io = x.astype(int)
        jo = y.astype(int)
        i1 = io + 1
        j1 = jo + 1

        s1 = x - io
        t1 = y - jo
        so = 1 - s1
        to = 1 - t1

        return (so*(to*f[io, jo]+t1*f[io, j1])+s1*(to*f[i1, jo]+t1*f[i1, j1]))

    def diffuse_explicit_2d(self, f):
        return self.dt*self.cd_2nd_2d(f, self.NU, self.NU)

    def diffuse_implicit_2d(self, fo, f, diff_coeff):
        return ((fo[1:-1, 1:-1] + (diff_coeff * self.dt)/(self.dx**2 * self.dy**2)*(self.dy**2 * (f[2:, 1:-1]) + self.dx**2 * (f[1:-1, 2:] + f[1:-1, :2])))/(1 + (2*diff_coeff * self.dt)/(self.dx**2 * self.dy)))

    def apply_pressure_2dX(self, p, c):
        return self.dt*self.cd_1st_2dX(p, c)

    def apply_pressure_2dY(self, p, c):
        return self.dt*self.cd_1st_2dY(p, c)

    def apply_force_2d(self, g):
        return self.dt*g[1:-1, 1:-1]

    def calc_source_2d(self, u, v):
        return (self.cd_1st_2dX(u, self.rho/self.dt)+self.cd_1st_2dY(v, self.rho/self.dt))

    def relax_presurre_poisson_2d(self, p, src):
        p [1:-1, 1:-1] = ((self.dy**2 * (p[2:, 1:-1]+p[:-2,1:-1])+self.dx**2*(p[1:-1,1:-1]+p[1:-1, :2])-self.dx**2*self.dy*src[1:-1, 1:-1])/(2*(self.dx**2 + self.dy**2)))
        return p

    def transform_pressure_poisson_2d(self, p, src):
        srcTrans = np.fft.fft2(src[1:-1, 1:-1])
        kx, ky = np.meshgrid(np.fft.fftfreq(self.Nx-2, d=self.dx),np.fft.fftfreq(self.Ny-2, d=self.dy))
        denom = 1.0/(4-2*np.cos(2*np.pi*kx*self.dx) - 2*np.cos(2*np.pi*ky*self.dy))
        denom[0, 0] = 0
        p[1:-1, 1:-1] = np.real_if_close(np.fft.ifft2(-srcTrans*denom*self.dx*self.dy))
        return p

    def mag_curl_term_2dX(self, u, v, Bx, By):
        return self.dt*(
            By[1:-1, 1:-1]*(u[1:-1, 2:]-u[1:-1, :2])/self.dy
            + u[1:-1, 1:-1]*(By[1:-1, 2:] - By[1:-1, :-2])/self.dy
            - Bx[1:-1, 1:-1]*(v[1:-1, 2:] - v[1:-1, :-2])/self.dy
            + v[1:-1, 1:-1]*(Bx[1:-1, 2:]-Bx[1:-1, :-2])/self.dy)

    def mag_curl_term_2dY(self, u, v, Bx, By):
        return self.dt*(
                - By[1:-1, 1:-1]*(u[2:, 1:-1]-u[:-2, 1:-1])/self.dx
                + u[1:-1, 1:-1]*(By[2:, 1:-1]-By[:2, 1:-1])/self.dx
                - Bx[1:-1, 1:-1]*(Bx[2:, 1:-1] - Bx[:-2, 1:-1])/self.dx)


def mag_diffuse_2d(self, fx, fy):
    return np.array([self.dt*self.central_diff_2nd_2d(fx, self.ETA/self.MU, self.ETA/self.MU), self.dt*self.central_diff_2nd_2d(fy, self.ETA/self.MU, self.ETA/self.MU)])

def mag_diffuse_2d(self,Bx,By):
    const = (self.ETA/self.MU*self.dt)/(self.dx**2*self.dy**2)
    return np.array([Bx[1:-1, 1:-1]/((1+2*const*(self.dy**2+self.dx**2)))+const/((1+2*const*(self.dy**2+self.dx**2)))*(self.dy**2*(Bx[2:, 1:-1]+Bx[:-2, 1:-1])+self.dx**2*(Bx[1:-1, 2:] + Bx[1:-1, :-2])),By[1:-1, 1:-1]/((1+2*const*(self.dy**2+self.dx**2)))+const/((1+2*const*(self.dy**2+self.dx**2)))*(self.dy**2*(By[2:,1:-1]+By[:-2,1:-1])+self.dx**2*(By[1:-1,2:1]+By[1:-1,1:-2]))])









# import required modules
from __future__ import division
import sys
import numpy as np
import pandas as pd
from scipy.optimize import fmin
from scipy import stats
from math import ceil
#from plotly.offline import plot
#import plotly.graph_objs as go
import csv
import os
import constants
import meshconstants


#########################################################
#            FINITE DIFFERENCE METHOD                   #
#                                                       #
#  A program for one dimensional flow in open channel   #
#                                                       #
#########################################################

#def compute_sections(section_arr, time_step):
def compute_sections(section_arr
        , j_current
        , j_next
        , upstream_flow_current
        , upstream_flow_next
        , downstream_stage_current
        , downstream_stage_next):
    #Set the upstream boundary flow
    #This is set as a delta -- a correction factor -- to work in the
    #predictor sweep.
    #So we set the Delta Q predictor, which is set from the input
    #time series for the upstream-most point.
    section_arr[0].time_steps[j_current].delta_flow_predictor = upstream_flow_next - upstream_flow_current

    for i, section in enumerate(section_arr):
        section_j = section.time_steps[j_current]
        section_j.depth = downstream_stage_current if i == 0 else section_j.water_z - section.bottom_z
        # print(f'depth {section_j.depth}')
        # print(f'bw {section.bottom_width}')
        # print(f'current stage {downstream_stage_current}')
        # print(f'next stage {downstream_stage_next}')
        # section_j.flow_area = section.bottom_width * section_j.depth
        section_j.flow_area = section.get_area_depth(section_j.depth)
        #TODO: All of these
        #Calculations based on the geometry should be
        #in class methods

                #TODO: not compatible with generalized sections.
        section_j.ci1 = section.bottom_width * (section_j.depth ** 2.0) / 2.0
        #section_j.hy = section_j.flow_area / (2.0 * section_j.depth + section.bottom_width)
        section_j.hy = section_j.flow_area / section.get_wetted_perimeter_depth(section_j.depth)
        section_j.conveyance = section.manning_k * section_j.flow_area * section_j.hy ** (2.0/3.0)
        # print (f'c0 {section_j.c0}')
        # print (f'manning_k {section.manning_k}')
        # print (f'flow_area {section_j.flow_area}')
        # print (f'Rw {section_j.hy}')
        #C       write(*,*)d(i),flow_area(i),ci1(i),y(j,i),z(i)
        #TODO: Where Does this WRITE TO???
        if i == 0: #treat the downstream boundary as a special case
            continue
        else:
            section_US = section_arr[i-1]
            section_US_j = section_US.time_steps[j_current]
                #TODO: not compatible with generalized sections.
            section_j.ci2_ds = ((section_j.depth ** 2.0) * (section_US_j.depth ** 2.0)) \
                          * (section.bottom_width - section_US.bottom_width) \
                          / (section_US.dx_ds * 4.0)
            section_US.bed_slope_ds = (section_US.bottom_z - section.bottom_z) \
                                    / section_US.dx_ds
            section_US_j.friction_slope_ds = section_US.manning_n_ds \
                           * 0.5 * section_US_j.flow \
                           * abs(section_US_j.flow) / (section_US_j.conveyance ** 2.0) \
                           + section_US.manning_n_ds * 0.5 * section_j.flow \
                           * abs(section_j.flow) / (section_j.conveyance ** 2.0)
            section_US_j.as0_ds = (section_US_j.flow_area + section_j.flow_area) \
                           / 2.0 * (section_US.bed_slope_ds \
                                    - section_US_j.friction_slope_ds)
            #TODO: This constants.GRAVITY will cause problems because it is not tied to any
            # definition of units
            section_US_j.gs0_ds = constants.GRAVITY * (section_US.bed_slope_ds \
                                    - section_US_j.friction_slope_ds)
            section_US.dbdx_ds = (section.bottom_width - section_US.bottom_width) \
                                  / section_US.dx_ds
    '''
      C in this version, we are computing the section attributes as they were
      C in the previous time step -- essentially as the initial condition for
      C the predictor step to advance to the corrector step and then to the
      C next time step
            subroutine section(n)
            parameter(grav=9.81)

            do 10 i=1,ncomp
            d(i)=y(n,i)-z(i)
            area(i)=bo(i)*d(i)
            ci1(i)=bo(i)*d(i)*d(i)/2.0
            hy=area(i)/(2.0*d(i)+bo(i))
            co(i)=sk(i)*area(i)*hy**(2.0/3.0)
      C       write(*,*)d(i),area(i),ci1(i),y(n,i),z(i)
      10    continue
      C       stop
            do 20 i=2,ncomp
            if(ityp(i-1).eq.1) then
            ci2(i)=(d(i)*d(i)+d(i-1)*d(i-1))*(bo(i)-bo(i-1))/(dx(i-1)*4.)
            beds=(z(i-1)-z(i))/dx(i-1)
            fs=f*0.5*q(n,i-1)*abs(q(n,i-1))/(co(i-1)*co(i-1))+
          1   f*0.5*q(n,i)*abs(q(n,i))/(co(i)*co(i))
            aso(i)=(area(i)+area(i-1))/2.0*(beds-fs)
            gso(i)=grav*(beds-fs)
            dbdx(i)=(bo(i)-bo(i-1))/dx(i-1)
            endif
      20    continue
            return
            end

    '''
    pass

def matrixp(section_arr
        , j_current
        , j_next
        , upstream_flow_current
        , upstream_flow_next
        , downstream_stage_current
        , downstream_stage_next):

    crmax = 0.0
    crmin = 100.0
    #TODO: crmax and crmin are not being used outside of this loop: make them useful
    for i, section in enumerate(section_arr):
        section_US = section_arr[i-1]
        #do 10 i=1,n_sections
        section_j = section.time_steps[j_current]
        #section_jnext = section.time_steps[j_next]
        section_j.velocity = section_j.flow / section_j.flow_area
        #u(i)=q(j,i)/flow_area(i)
        section_j.celerity = (constants.GRAVITY
                            * section_j.flow_area
                            / section.bottom_width) ** 0.5
        #c(i)=sqrt(grav*flow_area(i)/bo(i))

        #c     This is the matrix L (left eigenvector matrix - eq 13)
        e11 = 1.0
        if section_j.velocity == section_j.celerity: section_j.celerity = section_j.celerity + constants.CELERITY_EPSILON
        e12 = -1.0 / (section_j.velocity - section_j.celerity)
        e21 = 1.0
        e22 = -1.0 / (section_j.velocity + section_j.celerity)

        #c       L^{-1} (inverse of Left eigenvector matrix)
        f11 = -(section_j.velocity - section_j.celerity) / (2.0 * section_j.celerity)
        f12 = (section_j.velocity + section_j.celerity) / (2.0 * section_j.celerity)
        f21 = -(section_j.velocity ** 2.0 - section_j.celerity ** 2.0) / (2.0 * section_j.celerity)
        f22 = (section_j.velocity ** 2.0 - section_j.celerity ** 2.0) / (2.0 * section_j.celerity)

        #c       Diagonal wave matrix D (eq 12)
        d11 = abs(section_j.velocity + section_j.celerity)
        d22 = abs(section_j.velocity - section_j.celerity)
        #c       Equation 11 (L^{-1} D L)
        a11 = e11 * f11 * d11 + e21 * f12 * d22
        a12 = e12 * f11 * d11 + e22 * f12 * d22
        a21 = e11 * f21 * d11 + e21 * f22 * d22
        a22 = e12 * f21 * d11 + e22 * f22 * d22

        #dt(i)=dtini
        dt = 10
        #TODO: Make timestep handling more flexible by referring to jnext
        #dt = section_jnext.time - section_j.time

        #ter=bo(i)+2.0*flow_area(i)/bo(i)
        wetted_perimeter = section.get_wetted_perimeter_depth(section_j.depth)
        #dkda=sk(i)*((5.0/3.0*flow_area(i)**(2.0/3.0)*ter)-
        #1     (flow_area(i)**(5.0/3.0)*2.0/bo(i)))/ter/ter
        dkda = 1 / section.manning_n_ds \
            * (( 5.0 / 3.0 * section_j.flow_area ** (2.0/3.0) * wetted_perimeter) \
            - (section_j.flow_area ** (5.0/3.0) * 2.0 / section.bottom_width)) \
            / (wetted_perimeter ** 2.0)
        #TODO: not compatible with generalized sections.
        st11 = 0.0
        st12 = 0.0
        st21 = constants.GRAVITY * section_j.flow_area/section.bottom_width \
            / section.bottom_width * section.dbdx_ds + section_j.gs0_ds \
            + constants.MANNING_M * 2.0 * constants.GRAVITY \
            * section_j.flow_area * section_j.flow \
            * abs (section_j.flow) / section_j.conveyance ** 3.0 \
            * dkda
        st22 = -2.0 * constants.MANNING_M * section_j.flow_area \
             * section_j.flow_area \
             / section_j.conveyance

        if section.dx_ds == 0.0:
            section_j.cour = dt
        else:
            cour = dt / section.dx_ds
            crmax = max(crmax , cour * max(d11, d22))
            crmin = min(crmin , cour * max(d11, d22))

        # c     LHS of eq 7
        section_j.b11 = 0.5 - meshconstants.PHI - meshconstants.THETA  * cour * a11 - 0.5 * meshconstants.THETAS * st11 * dt
        section_j.b12 = -meshconstants.THETA  * cour * a12 - 0.5 * meshconstants.THETAS * st12 * dt
        section_j.b21 = -meshconstants.THETA * cour * a21 - 0.5 * meshconstants.THETAS * st21 * dt
        section_j.b22 = 0.5 - meshconstants.PHI - meshconstants.THETA  * cour * a22 - 0.5 * meshconstants.THETAS * st22 * dt
        # c
        if i == 1:
            cour = dt
        elif section_US.dx_ds == 0.0:
            cour = dt
        else:
            cour = dt / section_US.dx_ds
        # c
        g11 = 0.5 + meshconstants.PHI + meshconstants.THETA * cour * a11 - 0.5 * meshconstants.THETAS * st11 * dt
        g12 = meshconstants.THETA * cour * a12 - 0.5 * meshconstants.THETAS * st12 * dt
        g21 = meshconstants.THETA * cour * a21 - 0.5 * meshconstants.THETAS * st21 * dt
        g22 = 0.5 + meshconstants.PHI + meshconstants.THETA * cour * a22 - 0.5 * meshconstants.THETAS * st22 * dt
        # c
        section_j.g11inv =  g22 / (g11 * g22 - g12 * g21)
        section_j.g12inv = -g12 / (g11 * g22 - g12 * g21)
        section_j.g21inv = -g21 / (g11 * g22 - g12 * g21)
        section_j.g22inv =  g11 / (g11 * g22 - g12 * g21)
        #c
        section_j.f1 = section_j.flow
        section_j.f2 = section_j.flow ** 2.0 \
                    / section_j.flow_area \
                    + constants.GRAVITY * section_j.ci1
        #c
        # if not at the boundaries
        '''
        if (i >= 2.and.i.lt.I_UPSTREAM) then
            dip1=flow_area(i+1)/bo(i+1)
            di=2*flow_area(i)/bo(i)
            dim1=flow_area(i-1)/bo(i-1)
            eps2(i)=alfa2*abs(dip1-di+dim1)/(dip1+di+dim1)
        #c      write(*,*)i,eps2(i),flow_area(i),flow_area(i+1),flow_area(i-1)
        #endif
      #10    continue
      #c
    eps2(1) = eps2(2)
    eps2(I_UPSTREAM) = eps2(I_UPSTREAM-1)
      c
            do 20 i=2,n_sections-1
              if(ityp(i).ne.1) then
                eps2(i)=eps2(i-1)
                eps2(i+1)=eps2(i+2)
              endif
      20    continue
            do 40 i=2,n_sections-1
            eps2(i)=max(eps2(i+1),eps2(i))
      c      u(i)=(u(i+1)+u(i))/2.0
      c      c(i)=(c(i+1)+c(i))/2.0
            eps4(i)=max(0.,alfa4-eps2(i)/(u(i)+c(i)))
      c      write(*,*)i,eps2(i)
      40    continue
            d1(1)=0.0
            d2(1)=0.0
            d1(I_UPSTREAM)=0.0
            d2(I_UPSTREAM)=0.0
      c
            do 50 i=2,n_sections-1
            d11=abs(u(i)+c(i))
            d22=abs(u(i)-c(i))
            ei=max(d11,d22)
            d11=abs(u(i+1)+c(i+1))
            d22=abs(u(i+1)-c(i+1))
            ei1=max(d11,d22)
            eia=(ei+ei1)/2.0
            if(ityp(i).ne.1) then
              d1(i)=0.0
              d2(i)=0.0
            elseif(i.eq.2.or.i.eq.(I_UPSTREAM-1)) then ##TODO: Figure out if this is being treated correctly -- Would expect 'n_sections-1' to be more meaningful index, but the translation is confusing
              d1(i)=eps2(i)*eia*(flow_area(i+1)-flow_area(i))
              d2(i)=eps2(i)*eia*(q(j,i+1)-q(j,i))
      c      write(*,*)i,d1(i),d2(i),eps2(i),flow_area(i+1),flow_area(i)
            else
            d1(i)=eps2(i)*eia*(flow_area(i+1)-flow_area(i))-eps4(i)*(flow_area(i+2)-
          1       3*flow_area(i+1)+3*flow_area(i)-flow_area(i-1))
            d2(i)=eps2(i)*eia*(q(j,i+1)-q(j,i))-eps4(i)*(q(j,i+2)-
          1       3*q(j,i+1)+3*q(j,i)-q(j,i-1))
            endif
      c      write(*,*)i,d1(i),d2(i),eps2(i),flow_area(i+1),flow_area(i)
      50    continue
            return
            end

            '''

def apply_predictor(section_arr
        , j_current
        , j_next
        , upstream_flow_current
        , upstream_flow_next
        , downstream_stage_current
        , downstream_stage_next):
    pass

def dsbc():
    '''
      c
      c----------------------------------------------------------------------
      c     Downstream boundary condition
            subroutine dsbc(n)
      c----------------------------------------------------------------------
      c
            parameter(grav=9.81)
            common/arrays/ flow_area(1000),y(100000,1000),q(100000,1000),bo(1000),
          1    areap(1000),qp(1000),z(1000),delta_flow_predictor(1000),av11(1000),av12(1000),
          1                                       av21(1000),av22(1000),
          2              delta_flow_corrector(1000),delta_area_predictor(1000),delta_area_corrector(1000),ci1(1000),ci2(1000),
          3           aS0(1000),d(1000),f1(1000),g11inv(1000),g12inv(1000),
          4g21inv(1000),g22inv(1000),f2(1000),b11(1000),b12(1000),b21(1000),
          5b22(1000),eps2(1000),eps4(1000),d1(1000),d2(1000),u(1000),c(1000),
          6                        sk(1000),co(1000),gS0(1000),dbdx(1000),
          7                      dt(1000),ityp(1000),dx(1000)
            common/var/n_sections,cfl,f,qus,yds,rhs1,rhs2,c11,c12,c21,c22,us,time_step_optimization,
          1           thes,vv,dtini
            common/matrix/phi,theta,alfa2,alfa4,w
            common/sgate/ag,bg,option,yn,eps,dmeu,yw,bw,ydsn
      c
      c compute conjugate depth at dwon stream end
            ads=(y(j,I_UPSTREAM)-z(I_UPSTREAM))*bo(I_UPSTREAM)
            frds=q(j,I_UPSTREAM)/sqrt(grav*ads**3.0/bo(I_UPSTREAM))
            yconj=0.5*(y(j,I_UPSTREAM)-z(I_UPSTREAM))*(sqrt(1.0+8.0*frds*frds)-1.0)
            yconj=yconj+z(I_UPSTREAM)
            write(*,*)yn,y(j,I_UPSTREAM),yconj,z(I_UPSTREAM)
            if(yconj.lt.yn) then
      c      write(*,*)'no'
            if(option.eq.1.0) then
      c
      c downstream water level imposed (option 1)
                delta_area_corrector(I_UPSTREAM)=(yn-y(j,I_UPSTREAM))*bo(I_UPSTREAM)
                delta_area_predictor(I_UPSTREAM)=(yn-y(j,I_UPSTREAM))*bo(I_UPSTREAM)
                delta_flow_corrector(I_UPSTREAM)=delta_flow_predictor(I_UPSTREAM)
      c          write(*,*)'yes',delta_area_corrector(I_UPSTREAM),yn,y(j,I_UPSTREAM),delta_flow_corrector(I_UPSTREAM)
              elseif(option.eq.2.0) then
      c
      c downstream flow imposed  (option 2)
                delta_flow_corrector(I_UPSTREAM)=0.0
                delta_flow_predictor(I_UPSTREAM)=0.0
                delta_area_corrector(I_UPSTREAM)=delta_area_predictor(I_UPSTREAM)
              elseif(option.eq.3.0) then
      c
      c downstream rating curve imposed (option 3)
                delta_area_corrector(I_UPSTREAM)=delta_area_predictor(I_UPSTREAM)
                yn=(flow_area(I_UPSTREAM)+delta_area_predictor(I_UPSTREAM))/bo(I_UPSTREAM)
                qn=0.65*10*1.0*sqrt(2.0*grav*(yn-0.5))
                delta_flow_predictor(I_UPSTREAM)=qn-q(j,I_UPSTREAM)
                delta_flow_corrector(I_UPSTREAM)=delta_flow_predictor(I_UPSTREAM)
            endif
            else
      c
      c super critical flow exist at downstream exit
              delta_area_corrector(I_UPSTREAM)=delta_area_predictor(I_UPSTREAM)
              delta_flow_corrector(I_UPSTREAM)=delta_flow_predictor(I_UPSTREAM)
            endif
            return
            end
      c----------------------------------------------------------------------
    '''
    pass

def secpred(section_arr
        , j_current
        , j_next
        , upstream_flow_current
        , upstream_flow_next
        , downstream_stage_current
        , downstream_stage_next):
    '''
            subroutine secpred
            parameter(grav=9.81)
            common/arrays/ flow_area(1000),y(100000,1000),q(100000,1000),bo(1000),
          1    areap(1000),qp(1000),z(1000),delta_flow_predictor(1000),av11(1000),av12(1000),
          1                                       av21(1000),av22(1000),
          2              delta_flow_corrector(1000),delta_area_predictor(1000),delta_area_corrector(1000),ci1(1000),ci2(1000),
          3           aS0(1000),d(1000),f1(1000),g11inv(1000),g12inv(1000),
          4g21inv(1000),g22inv(1000),f2(1000),b11(1000),b12(1000),b21(1000),
          5b22(1000),eps2(1000),eps4(1000),d1(1000),d2(1000),u(1000),c(1000),
          6                        sk(1000),co(1000),gS0(1000),dbdx(1000),
          7                      dt(1000),ityp(1000),dx(1000)
            common/var/n_sections,cfl,f,qus,yds,rhs1,rhs2,c11,c12,c21,c22,us,time_step_optimization,
          1           thes,vv,dtini
            common/matrix/phi,theta,alfa2,alfa4,w
            common/sgate/ag,bg,option,yn,eps,dmeu,yw,bw,ydsn
      c
      c
            do 10 i=1,n_sections
            d(i)=areap(i)/bo(i)
            ci1(i)=bo(i)*d(i)*d(i)/2.0
            hy=areap(i)/(2.0*d(i)+bo(i))
            co(i)=sk(i)*areap(i)*hy**(2.0/3.0)
      10    continue
            do 20 i=1,n_sections-1
            if(ityp(i).eq.1) then
            ci2(i)=(d(i)*d(i)+d(i+1)*d(i+1))*(bo(i+1)-bo(i))/(dx(i)*4.)
            beds=(z(i)-z(i+1))/dx(i)
            fs=f*0.5*qp(i)*abs(qp(i))/(co(i)*co(i))+
          1   f*0.5*qp(i+1)*abs(qp(i+1))/(co(i+1)*co(i+1))
            aS0(i)=(areap(i)+areap(i+1))/2.0*(beds-fs)
            gS0(i)=grav*(beds-fs)
            dbdx(i)=(bo(i+1)-bo(i))/dx(i)
            endif
      20    continue
            return
            end
    '''
    pass


def matrixc(section_arr
        , j_current
        , j_next
        , upstream_flow_current
        , upstream_flow_next
        , downstream_stage_current
        , downstream_stage_next):
    '''
            subroutine matrixc
            parameter(grav=9.81)
            common/arrays/ flow_area(1000),y(100000,1000),q(100000,1000),bo(1000),
          1    areap(1000),qp(1000),z(1000),delta_flow_predictor(1000),av11(1000),av12(1000),
          1                                       av21(1000),av22(1000),
          2              delta_flow_corrector(1000),delta_area_predictor(1000),delta_area_corrector(1000),ci1(1000),ci2(1000),
          3           aS0(1000),d(1000),f1(1000),g11inv(1000),g12inv(1000),
          4g21inv(1000),g22inv(1000),f2(1000),b11(1000),b12(1000),b21(1000),
          5b22(1000),eps2(1000),eps4(1000),d1(1000),d2(1000),u(1000),c(1000),
          6                        sk(1000),co(1000),gS0(1000),dbdx(1000),
          7                      dt(1000),ityp(1000),dx(1000)
            common/var/n_sections,cfl,f,qus,yds,rhs1,rhs2,c11,c12,c21,c22,us,time_step_optimization,
          1           thes,vv,dtini
            common/matrix/phi,theta,alfa2,alfa4,w
            common/sgate/ag,bg,option,yn,eps,dmeu,yw,bw,ydsn
            c
            do 10 i=1,n_sections
            u(i)=qp(i)/areap(i)
            c(i)=sqrt(grav*areap(i)/bo(i))
      c
            e11=1.0
            if(u(i).eq.c(i)) c(i)=c(i)+0.00001
            e12=-1.0/(u(i)-c(i))
            e21=1.0
            e22=-1.0/(u(i)+c(i))
      c
            f11=-(u(i)-c(i))/(2.0*c(i))
            f12=(u(i)+c(i))/(2.0*c(i))
            f21=-(u(i)*u(i)-c(i)*c(i))/(2.0*c(i))
            f22=(u(i)*u(i)-c(i)*c(i))/(2.0*c(i))
      c
            d11=abs(u(i)+c(i))
            d22=abs(u(i)-c(i))
      c
            a11=e11*f11*d11+e21*f12*d22
            a12=e12*f11*d11+e22*f12*d22
            a21=e11*f21*d11+e21*f22*d22
            a22=e12*f21*d11+e22*f22*d22
      c
            dt(i)=dtini
      c
      c     Calculating dK/dA (eq 15)
            ter=bo(i)+2.0*areap(i)/bo(i)
            dkda=sk(i)*((5.0/3.0*areap(i)**(2.0/3.0)*ter)-
          1     (areap(i)**(5.0/3.0)*2.0/bo(i)))/ter/ter
      c
      c     Matrix S (eq 14)
            st11=0.0
            st12=0.0
            st21=grav*areap(i)/bo(i)/bo(i)*dbdx(i)+gS0(i)+f*
          1     2.0*grav*areap(i)*qp(i)*abs(qp(i))/co(i)**3.0*dkda
            st22=-2*f*qp(i)*areap(i)/co(i)/co(i)
      c
      c     cour == sigma
            if(i == 1) then
              cour=dt(i)
            else if(dx(i-1) == 0.0) then
              cour=dt(i)
            else
              cour=dt(i)/dx(i-1)
            endif
            b11(i)=0.5-phi-theta*cour*a11+0.5*thes*st11*dt(i)
            b12(i)=-theta*cour*a12+0.5*thes*st12*dt(i)
            b21(i)=-theta*cour*a21+0.5*thes*st21*dt(i)
            b22(i)=0.5-phi-theta*cour*a22+0.5*thes*st22*dt(i)
      c
            if(dx(i).eq.0.0) then
              cour=dt(i)
            else
              cour=dt(i)/dx(i)
              crmax=max(crmax,cour*max(d11,d22))
              crmin=min(crmin,cour*max(d11,d22))
            endif
            g11=0.5+phi+theta*cour*a11+0.5*thes*st11
            g12=theta*cour*a12+0.5*thes*st12
            g21=theta*cour*a21+0.5*thes*st21
            g22=0.5+phi+theta*cour*a22+0.5*thes*st22
      c
            g11inv(i)= g22/(g11*g22-g12*g21)
            g12inv(i)=-g12/(g11*g22-g12*g21)
            g21inv(i)=-g21/(g11*g22-g12*g21)
            g22inv(i)= g11/(g11*g22-g12*g21)
      c
            f1(i)=qp(i)
            f2(i)=qp(i)*qp(i)/areap(i)+grav*ci1(i)
      c
            if(i.ge.2.and.i.lt.n_sections) then
            dip1=areap(i+1)/bo(i+1)
            di=2*areap(i)/bo(i)
            dim1=areap(i-1)/bo(i-1)
            eps2(i)=alfa2*abs(dip1-di+dim1)/(dip1+di+dim1)
            endif
      10    continue
            eps2(1)=eps2(2)
            eps2(I_UPSTREAM)=eps2(I_UPSTREAM-1)
      c
            do 20 i=2,n_sections-1
              if(ityp(i).ne.1) then
                eps2(i)=eps2(i-1)
                eps2(i+1)=eps2(i+2)
              endif
      20    continue
      c
            do 40 i=n_sections-1,1,-1
            eps2(i+1)=max(eps2(i+1),eps2(i))
      c      u(i+1)=(u(i+1)+u(i))/2.0
      c      c(i+1)=(c(i+1)+c(i))/2.0
            eps4(i+1)=max(0.,alfa4-eps2(i+1)/(u(i+1)+c(i+1)))
      c      write(*,*)'corr',i,eps2(i)
      40    continue
            d1(1)=0.0
            d2(1)=0.0
            d1(I_UPSTREAM)=0.0
            d2(I_UPSTREAM)=0.0
      c
            do 50 i=2,n_sections-1
            d11=abs(u(i)+c(i))
            d22=abs(u(i)-c(i))
            ei=max(d11,d22)
            d11=abs(u(i-1)+c(i-1))
            d22=abs(u(i-1)-c(i-1))
            ei1=max(d11,d22)
            eia=(ei+ei1)/2.0
            if(ityp(i-1).ne.1) then
              d1(i)=0.0
              d2(i)=0.0
            elseif(i.eq.2.or.i.eq.(n_sections-1)) then
              d1(i)=eps2(i)*eia*(areap(i)-areap(i-1))
              d2(i)=eps2(i)*eia*(qp(i)-qp(i-1))
            else
            d1(i)=eps2(i)*eia*(areap(i)-areap(i-1))-eps4(i)*
          1        (areap(i+1)-3*areap(i)+3*areap(i-1)-areap(i-2))
            d2(i)=eps2(i)*eia*(qp(i)-qp(i-1))-eps4(i)*(qp(i+1)-3*qp(i)+
          1      3*qp(i-1)-qp(i-2))
            endif
      c      write(*,*)'corr',i,d1(i),d2(i),eps2(i),flow_area(i+1),flow_area(i)
      50    continue
            return
            end
    '''
    pass

def main():

    # loop on time and begin time-based simulation
    # This is the beginning of the real 'Guts' of the code
    for j in range(ntim-1):
        section(j) # TODO: Make this function
'''

        # There is the possibility that the Theta for the Predictor and
        # Corrector steps could be numerically different -- but we don't need to
        # take advantage of that here, so theta and theta_s are set to be
        # equal.
        thes = thetas
        matrixp(j) # TODO: Make this function

        for i in range(1, n_sections):
            sigma = dt[i] / dx[i-1]
            rhs1 =- sigma * (f1[i] - f1[i-1] - d1[i] + d1[i-1])
            rhs2 =- sigma * (f2[i] - f2[i-1] - d2[i] + d2[i-1])\
                + dt[i] * grav * (ci2[i] + aS0[i])
            c11=g11inv[i] * b11[i-1] + g12inv[i] * b21[i-1]
            c12=g11inv[i] * b12[i-1] + g12inv[i] * b22[i-1]
            c21=g21inv[i] * b11[i-1] + g22inv[i] * b21[i-1]
            c22=g21inv[i] * b12[i-1] + g22inv[i] * b22[i-1]
            delta_area_predictor[i] = g11inv[i] * rhs1 + g12inv[i] * rhs2\
                - c11*delta_area_predictor[i-1]-c12*delta_flow_predictor[i-1]
            delta_flow_predictor[i] = g21inv[i] * rhs1 + g22inv[i] * rhs2\
                - c21 * delta_area_predictor[i-1] - c22 * delta_flow_predictor[i-1]

  #       Boundary conditions at downstream (right boundary) -- section[I_UPSTREAM]
          delta_area_corrector(I_UPSTREAM)=delta_area_predictor(I_UPSTREAM)
          yn=(flow_area(I_UPSTREAM)+delta_area_predictor(I_UPSTREAM))/bo(I_UPSTREAM)
          arean=yn*bo(I_UPSTREAM)
          perimn=2.0*yn+bo(I_UPSTREAM)
          hyrdn=arean/perimn
          s0ds=-((z(I_UPSTREAM)-z(I_UPSTREAM-1))/dx(I_UPSTREAM))
          qn=skk*arean*hyrdn**(2.0/3.0)*sqrt(s0ds)
          delta_flow_predictor(I_UPSTREAM)=qn-q(j,I_UPSTREAM)
          delta_flow_corrector(I_UPSTREAM)=delta_flow_predictor(I_UPSTREAM)
        elseif(option_dsbc.eq.2)then
  c      	write(*,*) option_dsbc
          delta_area_corrector(I_UPSTREAM)=delta_area_predictor(I_UPSTREAM)
          yn=(flow_area(I_UPSTREAM)+delta_area_predictor(I_UPSTREAM))/bo(I_UPSTREAM)
          areac=yn*bo(I_UPSTREAM)
          perimc=2.0*yn+bo(I_UPSTREAM)
          hyrdc=areac/perimc
          s0ds=-((z(I_UPSTREAM)-z(I_UPSTREAM-1))/dx(I_UPSTREAM))
  c       qn=skk*areac*hyrdc**(2.0/3.0)*sqrt(s0ds)
          qcrit=1.05*(((yn**3.0)*(bo(I_UPSTREAM)**2.0)*grav)**(1.0/2.0))
          write(*,*)qcrit
          delta_flow_predictor(I_UPSTREAM)=qcrit-q(j,I_UPSTREAM)
          delta_flow_corrector(I_UPSTREAM)=delta_flow_predictor(I_UPSTREAM)
        else
          delta_area_corrector(I_UPSTREAM)=0.0
          delta_area_predictor(I_UPSTREAM)=0.0
          delta_flow_corrector(I_UPSTREAM)=delta_flow_predictor(I_UPSTREAM)
        endif
        '''

'''
  c      write(*,*)yn,arean,perimn,hyrdn
  c      write(*,*)qn,q(j,n_sections)
  c      stop
  c
  c     Update via predictor
        do 30 i=1,n_sections
        areap(i)=flow_area(i)+delta_area_predictor(i)
        qp(i)=q(j,i)+delta_flow_predictor(i)
  30    continue
  c
        call secpred
  c
        thes=thesinv
        call matrixc
  c
        do 40 i=n_sections-1,1,-1
          sigma=dt(i)/dx(i)
          rhs1=-sigma*(f1(i+1)-f1(i)-d1(i+1)+d1(i))
          rhs2=-sigma*(f2(i+1)-f2(i)-d2(i+1)+d2(i))
      1       +dt(i)*grav*(ci2(i)+aS0(i))
          c11=g11inv(i)*b11(i+1)+g12inv(i)*b21(i+1)
          c12=g11inv(i)*b12(i+1)+g12inv(i)*b22(i+1)
          c21=g21inv(i)*b11(i+1)+g22inv(i)*b21(i+1)
          c22=g21inv(i)*b12(i+1)+g22inv(i)*b22(i+1)
          delta_area_corrector(i)=g11inv(i)*rhs1+g12inv(i)*rhs2
      1       -c11*delta_area_corrector(i+1)-c12*delta_flow_corrector(i+1)
          delta_flow_corrector(i)=g21inv(i)*rhs1+g22inv(i)*rhs2
      1       -c21*delta_area_corrector(i+1)-c22*delta_flow_corrector(i+1)
  40    continue
  c
  c     Upstream boundary condition
  c     Prescribed flow at the upstream
  c     Area correction is calculated
        delta_flow_corrector(1)=delta_flow_predictor(1)
  c     delta_area_corrector(1)=delta_area_predictor(1)
  c
        do i=1,n_sections
  c       Final update
          da=(delta_area_predictor(i)+delta_area_corrector(i))/2.0
          dq=(delta_flow_predictor(i)+delta_flow_corrector(i))/2.0
          areanew=da+flow_area(i)
          if(areanew <= 0.0) areanew=0.0001
          y(j+1,i)=areanew/bo(i)+z(i)
          q(j+1,i)=q(j,i)+dq
        end do

        t = t + dtini
        print "('- cycle',i6,'  terminated')", j
        write(8, *) t, ((y(j+1, i) - z(i)) * bo(i), i=1,n_sections)
        write(9, *) t, (q(j+1, i), i=1,n_sections)

  10    continue

        close(8)
        close(9)

        stop
        end
'''

def apply_corrector(section_arr
        , j_current
        , j_next
        , upstream_flow_current
        , upstream_flow_next
        , downstream_stage_current
        , downstream_stage_next):
    pass

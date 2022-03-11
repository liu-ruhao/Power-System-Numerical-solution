#电力系统数值解法
import matplotlib.pyplot as plt
import config as C
import math as m

pi=m.pi

class Numerical_solution():
    def __init__(self,stru_config):
        print("初始条件为：")
        for i,v in stru_config.items():
            print("%s=%s"%(i,v))
        self.stru_config=stru_config 
        self.f=stru_config['f']
        self.Tj=stru_config['Tj']
        self.Pt=stru_config['Pt']
        self.P2m=stru_config['P2m']
        self.P3m=stru_config['P3m']
        self.y0=stru_config['y0']
        self.tm=stru_config['tm']
        self.h=stru_config['h']
        self.D=stru_config['D']
        
        self.train()


        
    def train(self):
        #计算极限角
        yh=pi-m.asin(self.Pt/self.P3m)
        cos_yc=(self.Pt*(yh-self.y0*pi/180)+self.P3m*m.cos(yh)-self.P2m*m.cos(self.y0*pi/180))/(self.P3m-self.P2m)
        yc=m.acos(cos_yc)*180/pi
        yh=yh*180/pi
        print('计算极限切除角')
        print('最大角是%s'%yh)
        print("极限切除角是%s°"%(yc))
        print("============================================")
        w0=360*self.f*pi/180 #w0的弧度
        xt,yt,wt=[0],[self.y0],[1]
        s,u,v,t=1,1,1,0
        
        print("计算极限切除时间:")
        for i in range(int(self.tm/self.h)):
            if t==0: #第一个时间段
                y=self.y0*pi/180 #0.60266                
                w_c_b=w0-w0 #0                      
                a_c_b=w0*(self.Pt-self.P2m*m.sin(y))/self.Tj   #deta的二阶导a,角加速度变化
                #时间段末速度估计值 v_last
                y_l=y+w_c_b*self.h
                w_l=w0+a_c_b*self.h
                #时间段末速度变化 v_chang_last
                w_c_l=w_l-w0
                a_c_l=pi/180*360*self.f*(self.Pt-self.P2m*m.sin(y_l))/self.Tj
                #时间段平均变化率 v_chang
                w_c=0.5*(w_c_b+w_c_l)
                a_c=0.5*(a_c_b+a_c_l)
                #时间段末速度 v
                y=(y+w_c*self.h)*180/pi
                w=w0+a_c*self.h
                t+=self.h
                xt.append(t)
                yt.append(y)
                wt.append(w/(self.f*360))
                #print("t=%s时，w：deta=%s : %s"%(t,w*180/pi/360/self.f,y))
            elif s:#故障，切割前
                #时间段开始速度的变化 v_change_begin
                y=y*pi/180                   #角度化为弧度，y=deta
                w_c_b=w-w0                            #deta的一阶导w,角速度变化
                a_c_b=w0*(self.Pt-self.P2m*m.sin(y))/self.Tj    #deta的二阶导a,角加速度变化
                #时间段末速度估计值 v_last
                y_l=y+w_c_b*self.h
                w_l=w+a_c_b*self.h
                #时间段末速度变化 v_chang_last
                w_c_l=w_l-w0
                a_c_l=w0*(self.Pt-self.P2m*m.sin(y_l))/self.Tj
                #时间段平均变化率 v_chang
                w_c=0.5*(w_c_b+w_c_l)
                a_c=0.5*(a_c_b+a_c_l)
                #时间段末速度 v
                y=(y+w_c*self.h)*180/pi
                w=w+a_c*self.h
                t+=self.h
                xt.append(t)
                yt.append(y)
                wt.append(w/w0)
                #print("t=%s时，w：deta=%s : %s"%(t,w/(self.f*360),y))
                if y>=yc:
                    s=0
                    y=y-w_c*self.h*180/pi
                    w=w-a_c*self.h
                    detay=w_c*self.h*180/pi
                    error=detay*100/yc
                    t-=self.h
                    xt.pop()
                    yt.pop()
                    wt.pop()
                    print("步长=%s秒时,极限切除角%s对应的估计极限切除时间为%s秒,"%(self.h,yc,t))
                    print("此时切割功率角为%s°,比极限切割功率角差%s"%(y,error)+'%')
                    print('此时角加速度为%s°'%(w*180/pi)+'/s')
                    print("可以设置更小的步长取得更精确的切割时间")
                    print("================================")
            #elif v:
            else:
                y=y*pi/180                   #角度化为弧度，y=deta
                w_c_b=w-w0                      #deta的一阶导w,角速度变化
                a_c_b=w0*(self.Pt-self.P3m*m.sin(y))/self.Tj   #deta的二阶导a,角加速度变化
                #时间段末速度估计值 v_last
                y_l=y+w_c_b*self.h
                w_l=w+a_c_b*self.h
                #时间段末速度变化 v_chang_last
                w_c_l=w_l-w0
                a_c_l=w0*(self.Pt-self.P3m*m.sin(y_l))/self.Tj
                #时间段平均变化率 v_chang
                w_c=0.5*(w_c_b+w_c_l)
                a_c=0.5*(a_c_b+a_c_l)
                #时间段末速度 v
                y=(y+w_c*self.h)*180/pi
                w=w+a_c*self.h
                t+=self.h
                xt.append(t)
                yt.append(y)
                wt.append(w/w0)
                #print("t=%s时，w：deta=%s : %s"%(t,w/(self.f*360),y))
                # if y<=(yh-130):
                    # v=0
                    # y=y-w_c*self.h*180/pi
                    # w=w-a_c*self.h
                    # detay=w_c*self.h*180/pi
                    # error=detay*100/yc
                    # t-=self.h
                    # xt.pop()
                    # yt.pop()
                    # wt.pop()
                    # print("步长=%s秒时,%s秒达到%s,比平衡功率角%s,差%s"%(self.h,t,y,(yh-90),error)+'/s')
                    # print("此时角速度为%s %s"%((w*180/pi),a_c))
                    # print("可以设置更小的步长取得更精确的切割时间")
                    # print("================================")
        ym=max(yt)
        i=yt.index(ym)
        th=xt[i]
        wh=wt[i]
        print("t=%s时,故障后达到的最大角是%s,相比理论最大角%s相差%s"%(th,ym,yh,abs(ym-yh)/yh*100),'%')
        print('此时角速度为%s'%(wt[i]*w0*180/pi))
        plt.plot(xt,yt)
        plt.xlabel('t/s')
        plt.ylabel('y/°')
        plt.plot(xt,wt)
        plt.show()

if __name__ == '__main__':
    main=Numerical_solution(C.stru_config)

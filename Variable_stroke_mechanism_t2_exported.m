classdef Variable_stroke_mechanism_t2_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        DRAGLINKCheckBox              matlab.ui.control.CheckBox
        CRANKROCKERCheckBox           matlab.ui.control.CheckBox
        TYPEFORLOOP1Label             matlab.ui.control.Label
        thetaFssEditField             matlab.ui.control.NumericEditField
        thetaFssEditFieldLabel        matlab.ui.control.Label
        thetaFcrsEditField            matlab.ui.control.NumericEditField
        thetaFcrsEditFieldLabel       matlab.ui.control.Label
        thetaFcsEditField             matlab.ui.control.NumericEditField
        thetaFcsEditFieldLabel        matlab.ui.control.Label
        thetaFbsEditField             matlab.ui.control.NumericEditField
        thetaFbsEditFieldLabel        matlab.ui.control.Label
        thetaFasEditField             matlab.ui.control.NumericEditField
        thetaFasEditFieldLabel        matlab.ui.control.Label
        FssEditField                  matlab.ui.control.NumericEditField
        FssEditFieldLabel             matlab.ui.control.Label
        FcrsEditField                 matlab.ui.control.NumericEditField
        FcrsEditFieldLabel            matlab.ui.control.Label
        FcsEditField                  matlab.ui.control.NumericEditField
        FcsEditFieldLabel             matlab.ui.control.Label
        FbsEditField                  matlab.ui.control.NumericEditField
        FbsEditFieldLabel             matlab.ui.control.Label
        FasEditField                  matlab.ui.control.NumericEditField
        FasEditFieldLabel             matlab.ui.control.Label
        crsEditField                  matlab.ui.control.NumericEditField
        crsEditFieldLabel             matlab.ui.control.Label
        csEditField                   matlab.ui.control.NumericEditField
        csEditFieldLabel              matlab.ui.control.Label
        bsEditField                   matlab.ui.control.NumericEditField
        bsEditFieldLabel              matlab.ui.control.Label
        asEditField                   matlab.ui.control.NumericEditField
        asEditFieldLabel              matlab.ui.control.Label
        NetaccelerationDropDown       matlab.ui.control.DropDown
        NetaccelerationDropDownLabel  matlab.ui.control.Label
        NetvelocityDropDown           matlab.ui.control.DropDown
        NetvelocityDropDownLabel      matlab.ui.control.Label
        ANALYSISLabel                 matlab.ui.control.Label
        AngularvelocityofDrivercrankradsEditField  matlab.ui.control.NumericEditField
        AngularvelocityofDrivercrankradsEditFieldLabel_2  matlab.ui.control.Label
        InclinationoffixedlinkdegreesEditField  matlab.ui.control.NumericEditField
        InclinationoffixedlinkdegreesEditFieldLabel  matlab.ui.control.Label
        Initialanglemadebycrank1withhorizontaldegreesEditField  matlab.ui.control.NumericEditField
        Initialanglemadebycrank1withhorizontaldegreesEditFieldLabel  matlab.ui.control.Label
        DurationofanimationsecsEditField  matlab.ui.control.NumericEditField
        DurationofanimationsecsEditFieldLabel  matlab.ui.control.Label
        TimestepsecsEditField         matlab.ui.control.NumericEditField
        TimestepsecsEditFieldLabel    matlab.ui.control.Label
        AnimationCheckBox             matlab.ui.control.CheckBox
        CLOSEButton                   matlab.ui.control.Button
        RESETButton                   matlab.ui.control.Button
        SAVEASButton                  matlab.ui.control.Button
        CALCULATEButton               matlab.ui.control.Button
        StaticforceanalysisCheckBox   matlab.ui.control.CheckBox
        VelocityandaccelerationanalysisCheckBox  matlab.ui.control.CheckBox
        ANGLEVARIATIONLabel           matlab.ui.control.Label
        TIMEINTERVALLabel             matlab.ui.control.Label
        LENGTHmEditField_6            matlab.ui.control.NumericEditField
        LENGTHmEditField_6Label       matlab.ui.control.Label
        ELabel                        matlab.ui.control.Label
        LENGTHmEditField_5            matlab.ui.control.NumericEditField
        LENGTHmEditField_5Label       matlab.ui.control.Label
        crLabel                       matlab.ui.control.Label
        LENGTHmEditField_4            matlab.ui.control.NumericEditField
        LENGTHmEditField_4Label       matlab.ui.control.Label
        dLabel                        matlab.ui.control.Label
        LENGTHmEditField_3            matlab.ui.control.NumericEditField
        LENGTHmEditField_3Label       matlab.ui.control.Label
        cLabel                        matlab.ui.control.Label
        LENGTHmEditField_2            matlab.ui.control.NumericEditField
        LENGTHmEditField_2Label       matlab.ui.control.Label
        bLabel                        matlab.ui.control.Label
        LENGTHmEditField              matlab.ui.control.NumericEditField
        LENGTHmEditFieldLabel         matlab.ui.control.Label
        aLabel                        matlab.ui.control.Label
        INPUTLabel                    matlab.ui.control.Label
    end

     methods (Access = public)
        
        function Animation_t2(app)               
                
                l1 =app.LENGTHmEditField.Value;
                l2 =app.LENGTHmEditField_2.Value;
                l3 = app.LENGTHmEditField_3.Value;
                l4  = app.LENGTHmEditField_4.Value; 
                
                l5 =app.LENGTHmEditField_5.Value; 
                offset = app.LENGTHmEditField_6.Value; 
                
                A = [l1;l2;l3;l4;l5];
                
                P=sort(A);
                
               if(P(1)+P(4)<P(2)+P(3))
                       
                       if(app.CRANKROCKERCheckBox.Value == 1)
                           a=P(1)
                           b=P(4);
                           c=P(2);
                           d=P(3);
                           cr=P(5);
                           e = offset;
                           s=sqrt((cr+c)^2-offset^2)-sqrt((cr-c)^2-offset^2);
                           if(cr-c+s>cr+c)
                               ch2 = 1;
                           else ch2 = 0;                           
                               
                           end
                       elseif(app.DRAGLINKCheckBox.Value ==1)
                           d=P(1);
                           a=P(4);
                           c=P(2);
                           b=P(3);
                           cr=P(5);
                           e = offset;
                           s=sqrt((cr+c)^2-offset^2)-sqrt((cr-c)^2-offset^2);
                           if(cr-c+s>cr+c)
                               ch2 = 1;
                           else ch2 = 0;
                           end
                           
                      end
                else ch2 = 0;                                        
                    
                end
                
                theta1=app.InclinationoffixedlinkdegreesEditField.Value;
                Wa=app.AngularvelocityofDrivercrankradsEditField.Value; %pi/10
                initial_angle=app.Initialanglemadebycrank1withhorizontaldegreesEditField.Value;
                time_of_running=app.DurationofanimationsecsEditField.Value;
                delt = app.TimestepsecsEditField.Value;
              
                final_angle=initial_angle+(time_of_running*Wa*180/pi);
                
                for theta=initial_angle:delt:final_angle
                
                k=((a^2-b^2+c^2+d^2)/2);
                A=-a.*(d.*cosd(theta1)-c).*cosd(theta)+k-(c*d*cosd(theta1))-a*d.*sind(theta).*sind(theta1);
                B=-(2.*a.*c.*sind(theta)-2*c*d.*sind(theta1));
                C=-a*(d*cosd(theta1)+c).*cosd(theta)+k+(c*d.*cosd(theta1))-a.*d.*sind(theta).*sind(theta1);
                si=2.*atand((-B-sqrt((B.^2)-4.*A.*C))./(2*A));
                
                G=((a^2+b^2-c^2+d^2)/2);
                D=-a.*(d.*cosd(theta1)+b).*cosd(theta)+G+(d*cosd(theta1)*b)-a*d*sind(theta1).*sind(theta);
                E=2*a*b.*sind(theta)-2*b*d.*sind(theta1);
                F=-(a.*(d*cosd(theta1)-b).*cosd(theta))+G-(b*d*cosd(theta1))-a.*d.*sind(theta).*sind(theta1);
                beta=2.*atand((-E+sqrt((E.^2)-4.*D.*F))./(2*D));                
                
                theta3=asind((e-c.*sind(si))./(cr));                
                
                plot([-2*a a+2*d+cr+e],[e e],'go-','LineWidth',1);hold on;
                
                plot([0 a*cosd(theta)], [0 a*sind(theta)],'ko-','LineWidth',2);hold on;
                plot([d*cosd(theta1) d*cosd(theta1)+c*cosd(si)], [d*sind(theta1) d*sind(theta1)+c*sind(si)], 'bo-','LineWidth',2); hold on;
                plot([0 d*cosd(theta1)], [0 d*sind(theta1)], 'mo-','LineWidth',2); hold on;
                plot([a*cosd(theta) a*cosd(theta)+b*cosd(beta)], [a*sind(theta) a*sind(theta)+b*sind(beta)], 'ro-','LineWidth',2);hold on;
                rectangle('position',[d+c*cosd(si)+cr*cosd(theta3)-P(1)/2 e-P(1)/2 P(1) P(1)],'FaceColor','c');hold on
                plot([d*cosd(theta1)+c*cosd(si) d+c*cosd(si)+cr*cosd(theta3)], [d*sind(theta1)+c*sind(si) e], 'yo-','LineWidth',2);hold off;
                grid on
                axis([-2*a a+2*d+cr+e -2*a a+2*d+cr+e]);
                pbaspect([1 1 1]);
                pause(0.01);
                drawnow
                end
                
            
        end
        
        function Velocity_accleration_At2(app)
                l1 =app.LENGTHmEditField.Value;
                l2 =app.LENGTHmEditField_2.Value;
                l3 = app.LENGTHmEditField_3.Value;
                l4  = app.LENGTHmEditField_4.Value; 
                
                l5 =app.LENGTHmEditField_5.Value; 
                offset = app.LENGTHmEditField_6.Value; 
                
                A = [l1;l2;l3;l4;l5];
                
                P=sort(A);
                
                    if(P(1)+P(4)<P(2)+P(3))
                       
                       if(app.CRANKROCKERCheckBox.Value == 1)
                           a=P(1)
                           b=P(4);
                           c=P(2);
                           d=P(3);
                           cr=P(5);
                           e = offset;
                           s=sqrt((cr+c)^2-offset^2)-sqrt((cr-c)^2-offset^2);
                           if(cr-c+s>cr+c)
                               ch2 = 1;
                           else ch2 = 0;                           
                               
                           end
                       elseif(app.DRAGLINKCheckBox.Value == 1)
                           d=P(1);
                           a=P(4);
                           c=P(2);
                           b=P(3);
                           cr=P(5);
                           e = offset;
                           s=sqrt((cr+c)^2-offset^2)-sqrt((cr-c)^2-offset^2);
                           if(cr-c+s>cr+c)
                               ch2 = 1;
                           else ch2 = 0;
                           end
                           
                      end
                    else ch2 = 0;                                        
                    
                    end
     
            theta1=app.InclinationoffixedlinkdegreesEditField.Value;
            Wa=app.AngularvelocityofDrivercrankradsEditField.Value; %pi/10
            initial_angle=app.Initialanglemadebycrank1withhorizontaldegreesEditField.Value;
            time_of_running=app.DurationofanimationsecsEditField.Value;
            delt = app.TimestepsecsEditField.Value;
            final_angle=initial_angle+(time_of_running*Wa*180/pi);
            
            theta=initial_angle:delt:final_angle;
            Aa = 0;
            time=linspace(0,time_of_running,length(theta));
            k=((a^2-b^2+c^2+d^2)/2);
            A=-a.*(d.*cosd(theta1)-c).*cosd(theta)+k-(c*d*cosd(theta1))-a*d.*sind(theta).*sind(theta1);
            B=-(2.*a.*c.*sind(theta)-2*c*d.*sind(theta1));
            C=-a*(d*cosd(theta1)+c).*cosd(theta)+k+(c*d.*cosd(theta1))-a.*d.*sind(theta).*sind(theta1);
            si=2.*atand((-B-sqrt((B.^2)-4.*A.*C))./(2*A));
            
            G=((a^2+b^2-c^2+d^2)/2);
            D=-a.*(d.*cosd(theta1)+b).*cosd(theta)+G+(d*cosd(theta1)*b)-a*d*sind(theta1).*sind(theta);
            E=2*a*b.*sind(theta)-2*b*d.*sind(theta1);
            F=-(a.*(d*cosd(theta1)-b).*cosd(theta))+G-(b*d*cosd(theta1))-a.*d.*sind(theta).*sind(theta1);
            beta=2.*atand((-E+sqrt((E.^2)-4.*D.*F))./(2*D));
            theta3=asind((e-c.*sind(si))./(cr));
            
            Wc=(a.*Wa.*sind(beta-theta))./(c.*sind(beta-si));
            Wb=(-a*Wa.*sind(si-theta))./(b.*sind(si-beta));
            Ab=((a.*Aa.*sind(si-theta))-(a.*(Wa.^2).*cosd(si-theta))-(b.*(Wb.^2).*cosd(si-beta))-(c.*(Wc.^2)))./(b.*sind(beta-si));
            Ac=((a.*Aa.*sind(beta-theta))-(a.*(Wa.^2).*cosd(beta-theta))-(b.*(Wb.^2))+(c.*(Wc.^2).*cosd(beta-si)))./(c.*sind(beta-si));
            Wcr=-(c*Wc.*cosd(si))./(cr.*cosd(theta3));
            v=((c*Wc.*sind(theta3-si))./(cosd(theta3)));
            
                       
            a3=(c*(Wc.^2).*sind(si)+cr.*(Wcr.^2).*sind(theta3)-c.*Ac.*cos(theta3))./cr.*cosd(theta3);
            a_piston=(c.*Ac.*sind(theta3-si)-c.*(Wc.^2).*cosd(theta3-si)-cr.*(Wcr.^2))./cosd(theta3);           
            
            
            %positions of each point
            O=[0;0];
            A=[a*cosd(theta);a*sind(theta)];
            B=[d*cosd(theta1)+c*cosd(si);d*sind(theta1)+c*sind(si)];
            
            A_x=A(1,:);
            A_y=A(2,:);
            
            A_vx=diff(A_x)./diff(time);
            A_vy=diff(A_y)./diff(time);
            
            A_v=sqrt(A_vx.^2+A_vy.^2);
            A_v=[0 A_v];
            
            A_a=diff(A_v)./diff(time);
            A_a=[0 A_a];            
                        
            B_x=B(1,:);
            B_y=B(2,:);
            
            B_vx=diff(B_x)./diff(time);
            B_vy=diff(B_y)./diff(time);
            
            B_v=sqrt(B_vx.^2+B_vy.^2);
            B_v=[0 B_v];
            
            B_a=diff(B_v)./diff(time);
            B_a=[0 B_a];    
                                              
                   
           if app.NetvelocityDropDown.Value
                value3 = app.NetvelocityDropDown.Value;
                   if strcmpi(value3,'Point A')                         
                         plot(time,A_v);grid on
                         ylabel('Velocity of point A (m/s)')
                         xlabel('time (secs)')
                   elseif strcmpi(value3,'Point B')                        
                        plot(time,B_v);grid on
                         ylabel('Velocity of point B (m/s)')
                         xlabel('time (secs)')
                   elseif strcmpi(value3,'Slider')                        
                           plot(time,v)
                           grid on
                           xlabel('time (secs)')
                           ylabel('velocity of slider (m/s)')
                   end
           end
              
                   
           if app.NetaccelerationDropDown.Value
               value4 = app.NetaccelerationDropDown.Value;
                   if strcmpi(value4,'Point A')                         
                         plot(time,A_a); grid on
                         ylabel('acceleration of point A (m/s^2)')
                         xlabel('time (secs)')
                   elseif strcmpi(value4,'Point B')                       
                       plot(time,B_a); grid on
                       ylabel('acceleration of point B (m/s^2)')
                       xlabel('time (secs)')
                   elseif strcmpi(value4,'Slider')                       
                       plot(time,a_piston)
                        grid on
                        xlabel('time (secs)')
                        ylabel('acceleration of the slider (m/s^2)')
                         
                   end
           end
     end
        
        
        
        function [ch] = check_conditions_t2(app)
                l1 =app.LENGTHmEditField.Value;
                l2 =app.LENGTHmEditField_2.Value;
                l3 = app.LENGTHmEditField_3.Value;
                l4  = app.LENGTHmEditField_4.Value; 
                
                l5 =app.LENGTHmEditField_5.Value; 
                offset = app.LENGTHmEditField_6.Value; 
                
                A = [l1;l2;l3;l4;l5];
                
                P=sort(A);
                
                     if(P(1)+P(4)<P(2)+P(3))
                       
                       if(app.CRANKROCKERCheckBox.Value == 1)
                           a=P(1)
                           b=P(4);
                           c=P(2);
                           d=P(3);
                           cr=P(5);
                           e = offset;
                           s=sqrt((cr+c)^2-offset^2)-sqrt((cr-c)^2-offset^2);
                           if(cr-c+s>cr+c)
                               ch2 = 1;
                           else ch2 = 0;                           
                               
                           end
                       elseif(app.DRAGLINKCheckBox.Value == 1)
                           d=P(1);
                           a=P(4);
                           c=P(2);
                           b=P(3);
                           cr=P(5);
                           e = offset;
                           s=sqrt((cr+c)^2-offset^2)-sqrt((cr-c)^2-offset^2);
                           if(cr-c+s>cr+c)
                               ch2 = 1;
                           else ch2 = 0;
                           end
                           
                      end
                else ch2 = 0;
                      ch=0;
                      return;                    
                    
                end
        
            theta1=app.InclinationoffixedlinkdegreesEditField.Value;
            Wa=app.AngularvelocityofDrivercrankradsEditField.Value; %pi/10
            initial_angle=app.Initialanglemadebycrank1withhorizontaldegreesEditField.Value;
            time_of_running=app.DurationofanimationsecsEditField.Value;
            delt = app.TimestepsecsEditField.Value;
            
            final_angle=initial_angle+time_of_running*Wa;
            if(ch2==1)
             for theta=initial_angle:delt:final_angle;
                
                k=((a^2-b^2+c^2+d^2)/2);
                A=-a.*(d.*cosd(theta1)-c).*cosd(theta)+k-(c*d*cosd(theta1))-a*d.*sind(theta).*sind(theta1);
                B=-(2.*a.*c.*sind(theta)-2*c*d.*sind(theta1));
                C=-a*(d*cosd(theta1)+c).*cosd(theta)+k+(c*d.*cosd(theta1))-a.*d.*sind(theta).*sind(theta1);
                si=2.*atand((-B-sqrt((B.^2)-4.*A.*C))./(2*A));
                
                G=((a^2+b^2-c^2+d^2)/2);
                D=-a.*(d.*cosd(theta1)+b).*cosd(theta)+G+(d*cosd(theta1)*b)-a*d*sind(theta1).*sind(theta);
                E=2*a*b.*sind(theta)-2*b*d.*sind(theta1);
                F=-(a.*(d*cosd(theta1)-b).*cosd(theta))+G-(b*d*cosd(theta1))-a.*d.*sind(theta).*sind(theta1);
                beta=2.*atand((-E+sqrt((E.^2)-4.*D.*F))./(2*D));
                theta3=asind((e-c.*sind(si))./(cr));
             end
            end        
 
            if  ((imag(beta) || imag(si) || imag(theta3)) )
                ch1 = (~((imag(beta) || imag(si) || imag(theta3)) ));
                ch = ch1 && ch2;          
            else  ch1 = 1; 
                 ch = ch1 && ch2;
            
            end
            return;
        end
        
        
        
        function [res2]=Staticforce_At2(app)
            
                l1 =app.LENGTHmEditField.Value;
                l2 =app.LENGTHmEditField_2.Value;
                l3 = app.LENGTHmEditField_3.Value;
                l4  = app.LENGTHmEditField_4.Value; 
                
                l5 =app.LENGTHmEditField_5.Value; 
                offset = app.LENGTHmEditField_6.Value; 
                
                A = [l1;l2;l3;l4;l5];
                
                P=sort(A);
                
                    if(P(1)+P(4)<P(2)+P(3))
                       
                       if(app.CRANKROCKERCheckBox.Value == 1)
                           a=P(1)
                           b=P(4);
                           c=P(2);
                           d=P(3);
                           cr=P(5);
                           e = offset;
                           s=sqrt((cr+c)^2-offset^2)-sqrt((cr-c)^2-offset^2);
                           if(cr-c+s>cr+c)
                               ch2 = 1;
                           else ch2 = 0;                           
                               
                           end
                       elseif(app.DRAGLINKCheckBox.Value == 1)
                           d=P(1);
                           a=P(4);
                           c=P(2);
                           b=P(3);
                           cr=P(5);
                           e = offset;
                           s=sqrt((cr+c)^2-offset^2)-sqrt((cr-c)^2-offset^2);
                           if(cr-c+s>cr+c)
                               ch2 = 1;
                           else ch2 = 0;
                           end
                           
                      end
                   else ch2 = 0;                                        
                    
                   end
            
            theta1=app.InclinationoffixedlinkdegreesEditField.Value;
            Wa=app.AngularvelocityofDrivercrankradsEditField.Value; %pi/10
            initial_angle=app.Initialanglemadebycrank1withhorizontaldegreesEditField.Value;
            time_of_running=app.DurationofanimationsecsEditField.Value;
            delt = app.TimestepsecsEditField.Value;
            
            final_angle=initial_angle+time_of_running*Wa;
            
            as = app.asEditField.Value;
            bs = app.bsEditField.Value;
            cs = app.csEditField.Value;
            crs = app.crsEditField.Value;
            Fas = app.FasEditField.Value;
            Fbs = app.FbsEditField.Value;
            Fcrs = app.FcrsEditField.Value;
            Fcs = app.FcsEditField.Value;
            Fss = app.FssEditField.Value;
            thetaFcs=app.thetaFcsEditField.Value;
            thetaFbs=app.thetaFbsEditField.Value;
            thetaFcrs=app.thetaFcrsEditField.Value;
            thetaFss = app.thetaFssEditField.Value;
            thetaFas=app.thetaFasEditField.Value;
            
            if as>a
                 f1= errordlg("Analysis not possible");
                 z=0;
            elseif bs>b
                 f1= errordlg("Analysis not possible");
                 z=0;
            elseif cs>c
                 f1= errordlg("Analysis not possible");
                 z=0;
            elseif crs>cr
                 f1= errordlg("Analysis not possible");
                 z=0;
            else  z=1;
                
            end
            
            theta=app.Initialanglemadebycrank1withhorizontaldegreesEditField.Value; %defined by user
            if(z==1)
            k=((a^2-b^2+c^2+d^2)/2);
            A=-a.*(d.*cosd(theta1)-c).*cosd(theta)+k-(c*d*cosd(theta1))-a*d.*sind(theta).*sind(theta1);
            B=-(2.*a.*c.*sind(theta)-2*c*d.*sind(theta1));
            C=-a*(d*cosd(theta1)+c).*cosd(theta)+k+(c*d.*cosd(theta1))-a.*d.*sind(theta).*sind(theta1);
            si=2.*atand((-B-sqrt((B.^2)-4.*A.*C))./(2*A));
            
            G=((a^2+b^2-c^2+d^2)/2);
            D=-a.*(d.*cosd(theta1)+b).*cosd(theta)+G+(d*cosd(theta1)*b)-a*d*sind(theta1).*sind(theta);
            E=2*a*b.*sind(theta)-2*b*d.*sind(theta1);
            F=-(a.*(d*cosd(theta1)-b).*cosd(theta))+G-(b*d*cosd(theta1))-a.*d.*sind(theta).*sind(theta1);
            beta=2.*atand((-E+sqrt((E.^2)-4.*D.*F))./(2*D));
            
            
            theta3=asind((e-c.*sind(si))./(cr));
            
            syms Fab 
            moment_1=cross(as*[cosd(theta) sind(theta) 0],Fas*[cosd(thetaFas) sind(thetaFas) 0])+cross(a*[cosd(theta) sind(theta) 0],Fab*[cosd(beta) sind(beta) 0])==0
            moment_1 = moment_1(3);
            Fab = double(rhs((isolate(moment_1,Fab))))
            Torque_1 = cross( a*[cosd(theta) sind(theta) 0], Fab*[cosd(beta) sind(beta) 0])
            fprintf("%f Nm ",Torque_1)
            
            
            syms Fab
            moment_2=cross([bs*cosd(beta) bs*sind(beta) 0],Fbs*[cosd(thetaFbs) sind(thetaFbs) 0])+cross([b*cosd(beta) b*sind(beta) 0],Fab*[cosd(si) sind(si) 0])==0
            moment_2 = moment_2(3);
            Fab = double(rhs((isolate(moment_2,Fab))))
            Torque_2 = cross( a*[cosd(theta) sind(theta) 0], Fab*([cosd(beta) sind(beta) 0]))
            fprintf("%f Nm ",Torque_2)
            
            syms Fab 
            moment_3=cross(cs*[cosd(si) sind(si) 0],Fcs*[cosd(thetaFcs) sind(thetaFcs) 0])+cross([c*cosd(si) c*sind(si) 0],Fab*[cosd(beta) sind(beta) 0])==0
            moment_3 = moment_3(3);
            Fab = double(rhs((isolate(moment_3,Fab))))
            Torque_3 = cross( a*[cosd(theta) sind(theta) 0], Fab*([cosd(si) sind(si) 0]))
            fprintf("%f Nm ",Torque_3)
            
            syms Fab  F
            Fx = Fss + F*cosd(theta3) ==0;
            Fbd = double(rhs(vpa(isolate(Fx,F))))
            moment_4=cross(c*[cosd(si) sind(si) 0],Fbd*[cosd(theta3) sind(theta3) 0])+cross(c*[cosd(si) sind(si) 0],Fab*[cosd(beta) sind(beta) 0])==0
            moment_4 = moment_4(3);
            Fab = double(rhs((isolate(moment_4,Fab))))
            Torque_4 = cross( a*[cosd(theta) sind(theta) 0], Fab*([cosd(beta) sind(beta) 0]))
            fprintf("%f Nm ",Torque_4)
            
            NET_TORQUE=Torque_1+Torque_2+Torque_3+Torque_4
            fprintf("%f Nm ",NET_TORQUE)
            res2 = NET_TORQUE(3);
            end
            return;
        end
        
        
        function create_outputfile(app)
                l1 =app.LENGTHmEditField.Value;
                l2 =app.LENGTHmEditField_2.Value;
                l3 = app.LENGTHmEditField_3.Value;
                l4  = app.LENGTHmEditField_4.Value; 
                
                l5 =app.LENGTHmEditField_5.Value; 
                offset = app.LENGTHmEditField_6.Value; 
                
                A = [l1;l2;l3;l4;l5];
                
                
                P=sort(A);
                
                      if(P(1)+P(4)<P(2)+P(3))
                       
                       if(app.CRANKROCKERCheckBox.Value == 1)
                           a=P(1)
                           t_1 = 'CRANK';
                           b=P(4);
                           t_2 = 'COUPLER';
                           c=P(2);
                           t_3 = 'ROCKER';
                           d=P(3);
                           t_4 = 'FIXED LINK';
                           cr=P(5);
                           t_5 = 'CONNECTING ROD';
                           e = offset;
                           t_6 = 'OFFSET';
                           s=sqrt((cr+c)^2-offset^2)-sqrt((cr-c)^2-offset^2);
                           if(cr-c+s>cr+c)
                               ch2 = 1;
                           else ch2 = 0;                           
                               
                           end
                       elseif(app.DRAGLINKCheckBox.Value == 1)
                           d=P(1);
                           t_1 = 'FIXED LINK'
                           a=P(4);
                           t_2 = 'CRANK'
                           c=P(2);
                           t_3 = 'ROCKER';
                           b=P(3);
                           t_4 = 'COUPLER';
                           cr=P(5);
                           t_5 = 'CONNECTING ROD';
                           e = offset;
                           t_6 = 'OFFSET';
                           s=sqrt((cr+c)^2-offset^2)-sqrt((cr-c)^2-offset^2);
                           if(cr-c+s>cr+c)
                               ch2 = 1;
                           else ch2 = 0;
                           end
                           
                      end
                     else ch2 = 0;                                        
                    
                     end
            
            theta1=app.InclinationoffixedlinkdegreesEditField.Value;
            Wa=app.AngularvelocityofDrivercrankradsEditField.Value; %pi/10
            initial_angle=app.Initialanglemadebycrank1withhorizontaldegreesEditField.Value;
            time_of_running=app.DurationofanimationsecsEditField.Value;
            delt = app.TimestepsecsEditField.Value;
            final_angle=initial_angle+(time_of_running*Wa*180/pi);
            
            theta=initial_angle:delt:final_angle;
            Aa = 0;
            time=linspace(0,time_of_running,length(theta));
            k=((a^2-b^2+c^2+d^2)/2);
            A=-a.*(d.*cosd(theta1)-c).*cosd(theta)+k-(c*d*cosd(theta1))-a*d.*sind(theta).*sind(theta1);
            B=-(2.*a.*c.*sind(theta)-2*c*d.*sind(theta1));
            C=-a*(d*cosd(theta1)+c).*cosd(theta)+k+(c*d.*cosd(theta1))-a.*d.*sind(theta).*sind(theta1);
            si=2.*atand((-B-sqrt((B.^2)-4.*A.*C))./(2*A));
            
            G=((a^2+b^2-c^2+d^2)/2);
            D=-a.*(d.*cosd(theta1)+b).*cosd(theta)+G+(d*cosd(theta1)*b)-a*d*sind(theta1).*sind(theta);
            E=2*a*b.*sind(theta)-2*b*d.*sind(theta1);
            F=-(a.*(d*cosd(theta1)-b).*cosd(theta))+G-(b*d*cosd(theta1))-a.*d.*sind(theta).*sind(theta1);
            beta=2.*atand((-E+sqrt((E.^2)-4.*D.*F))./(2*D));
            theta3=asind((e-c.*sind(si))./(cr));
            
            Wc=(a.*Wa.*sind(beta-theta))./(c.*sind(beta-si));
            Wb=(-a*Wa.*sind(si-theta))./(b.*sind(si-beta));
            Ab=((a.*Aa.*sind(si-theta))-(a.*(Wa.^2).*cosd(si-theta))-(b.*(Wb.^2).*cosd(si-beta))-(c.*(Wc.^2)))./(b.*sind(beta-si));
            Ac=((a.*Aa.*sind(beta-theta))-(a.*(Wa.^2).*cosd(beta-theta))-(b.*(Wb.^2))+(c.*(Wc.^2).*cosd(beta-si)))./(c.*sind(beta-si));
            Wcr=-(c*Wc.*cosd(si))./(cr.*cosd(theta3));
            v=((c*Wc.*sind(theta3-si))./(cosd(theta3)));
            
                       
            a3=(c*(Wc.^2).*sind(si)+cr.*(Wcr.^2).*sind(theta3)-c.*Ac.*cos(theta3))./cr.*cosd(theta3);
            a_piston=(c.*Ac.*sind(theta3-si)-c.*(Wc.^2).*cosd(theta3-si)-cr.*(Wcr.^2))./cosd(theta3);           
            
            
            %positions of each point
            O=[0;0];
            A=[a*cosd(theta);a*sind(theta)];
            B=[d*cosd(theta1)+c*cosd(si);d*sind(theta1)+c*sind(si)];
            
            A_x=A(1,:);
            A_y=A(2,:);
            
            A_vx=diff(A_x)./diff(time);
            A_vy=diff(A_y)./diff(time);
            
            A_v=sqrt(A_vx.^2+A_vy.^2);
            A_v=[0 A_v];
            
            A_a=diff(A_v)./diff(time);
            A_a=[0 A_a];            
                        
            B_x=B(1,:);
            B_y=B(2,:);
            
            B_vx=diff(B_x)./diff(time);
            B_vy=diff(B_y)./diff(time);
            
            B_v=sqrt(B_vx.^2+B_vy.^2);
            B_v=[0 B_v];
            
            B_a=diff(B_v)./diff(time);
            B_a=[0 B_a]; 
            
            as = app.asEditField.Value;
            bs = app.bsEditField.Value;
            cs = app.csEditField.Value;
            crs = app.crsEditField.Value;
            Fas = app.FasEditField.Value;
            Fbs = app.FbsEditField.Value;
            Fcrs = app.FcrsEditField.Value;
            Fcs = app.FcsEditField.Value;
            Fss = app.FssEditField.Value;
            thetaFcs=app.thetaFcsEditField.Value;
            thetaFbs=app.thetaFbsEditField.Value;
            thetaFcrs=app.thetaFcrsEditField.Value;
            thetaFss = app.thetaFssEditField.Value;
            thetaFas=app.thetaFasEditField.Value;

            
            tor2 = Staticforce_At2(app);
            
            startingFolder = userpath;
            filter = {'.docx';'*.dat';'*.txt';'*.m';'*.slx';'*.mat';'*.*'};
            defaultFileName = fullfile(startingFolder, filter);
            [baseFileName, folder] = uiputfile(defaultFileName, 'Specify a file');
            if baseFileName == 0
              % User clicked the Cancel button.
              return;
            end
            fullFileName = fullfile(folder, baseFileName);
            import mlreportgen.dom.*
            d_2 = Document(fullFileName,"docx");
            open(d_2);
            if d_2 ~= -1
                %1st given data table
               tableStyle = { Width("110%"), ...
                               Border("solid"), ...
                               RowSep("solid"), ...
                               ColSep("solid") };

               append(d_2,Heading1("GIVEN DATA:"));
                BodyContent = {'Time range (secs)', time_of_running; ...
                               'Time step (secs)', delt; ...
                               'Initial angle made by crank1 with horizontal (degrees)',initial_angle; ...
                                'Angle made by fixed link with horizontal (degrees)' ,theta1; ...                                 
                                'Anglular velocity of driver crank1 (rad/s) ' ,Wa};
                
                    tableContent_1 = [BodyContent];
                    
                    table = Table(tableContent_1);
                    table.Style = tableStyle;
                    
                    table.TableEntriesHAlign = "center";
                    append(d_2, table); 
                    
                    % 2nd - inputs table
                    tableStyle_1 = { Width("70%"), ...
                                   Border("solid"), ...
                                   RowSep("solid"), ...
                                   ColSep("solid") };
    
                    append(d_2,Heading1("INPUTS: "));
                    
                    HeaderContent = {'LINKS','LENGTH (metres)','TYPE'};
                    BodyContent = {'a',a,t_1;'b',b,t_2;'c',c,t_3;...
                                    'd',d,t_4;'cr',cr,t_5;'e',e,t_6;};
                    tableContent_2 = [HeaderContent;BodyContent];
                    
                    table = Table(tableContent_2);
                    table.Style = tableStyle_1;
                    
                    table.TableEntriesHAlign = "center";
                    append(d_2, table);
                    
                     %3rd inputs table
                    append(d_2,Heading1("FORCES APPLIED ON LINKS"));
                    
                    formatSpec = 'On CRANK %d (N) at %d (m) at an angle %d (deg) with horizontal';                   
                    str2 = sprintf(formatSpec,Fas,as,thetaFas);
                    append(d_2,str2);
                    
                    formatSpec = 'On COUPLER %d (N) at %d (m) at an angle %d (deg) with horizontal';                   
                    str2 = sprintf(formatSpec,Fbs,bs,thetaFbs);
                    append(d_2,str2);
                    
                    formatSpec = 'On ROCKER %d (N) at %d (m) at an angle %d (deg) with horizontal';                   
                    str2 = sprintf(formatSpec,Fcs,cs,thetaFcs);
                    append(d_2,str2);
                    
                                     
                    formatSpec = 'On CONNECTING ROD %d (N) at %d (m) at an angle %d (deg) with horizontal';                   
                    str2 = sprintf(formatSpec,Fcrs,crs,thetaFcrs);
                    append(d_2,str2);  
                    
                    formatSpec = 'On slider %d (N) at an angle %d (deg) with horizontal';                   
                    str2 = sprintf(formatSpec,Fss,thetaFss);
                    append(d_2,str2);                   
         
         
                    
                    %3rd - outputs table 1
                    
                    headerContent = {'Time (secs)','CRANK angle with horizontal (deg)',...
                        'ROCKER angle(deg)','COUPLER angle(deg)','COUPLER2 angle(deg)'};
                    bodyContent = [time',theta',si',beta',theta3'];
                    
                    data_str = string(bodyContent);
                    %round to 2 decimal places
                    for i = 1:numel(data_str)
                    data_str(i) = sprintf('%.2f',data_str(i));
                    end            
                    
                    tableContent = [headerContent; data_str];            
                    
                    append(d_2,Heading1("All Table Entries Centered"));
                    
                    table = Table(tableContent);
                    table.Style = tableStyle;
                    
                    table.TableEntriesHAlign = "center";
                    append(d_2, table);
            
            
                   % 3rd outputs table 2
                    table_Style = { Width("100%"), ...
                                   Border("solid"), ...
                                   RowSep("solid"), ...
                                   ColSep("solid") };
                    
                    header_Content = {'Time (secs)','Velocity_A (m/s)','Velocity_B (m/s)','Velocity of slider (m/s)',...
                              'Acceleration_A (m/s^2)','Acceleration_B (m/s^2)','Acceleration of slider (m/s^2)' };
                    body_Content = [time',A_v',B_v',v',A_a',B_a',a_piston'];
                    
                    data__str = string(body_Content);
                    %round to 2 decimal places
                    for i = 1:numel(data__str)
                    data__str(i) = sprintf('%.2f',data__str(i));
                    end
                    
                    
                    table_Content = [header_Content; data__str];
                    
                    
                    append(d_2,Heading1("All Table Entries Centered"));
                    
                    table = Table(table_Content);
                    table.Style = table_Style;
                    
                    table.TableEntriesHAlign = "center";
                    append(d_2, table);
                    
                    formatSpec = 'Net torque acting on crank1 is %d N/m';
                    A2 = tor2;
                    
                    str1 = sprintf(formatSpec,A2)
                    
                    append(d_2,Heading1(str1));
                    
                    close(d_2);                           
                
                
             else
                warningMessage = sprintf('Cannot open file:\n', fullFileName);
                uiwait(warndlg(warningMessage));           
  
        end
    end
 end  

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: CALCULATEButton
        function CALCULATEButtonPushed(app, event)
            val = check_conditions_t2(app);
            if val 
                if app.VelocityandaccelerationanalysisCheckBox.Value
                     Velocity_accleration_At2(app);
                elseif app.StaticforceanalysisCheckBox.Value
                     st2 = Staticforce_At2(app);
                     message = sprintf('Net torque acting on crank at equilibrium %d N/m',st2);
                     uiwait(msgbox(message,'modal'));
                 elseif app.AnimationCheckBox.Value
                     Animation_t2(app);
                end
            else f = errordlg('The given inputs do not the match the conditions required for the mechanism','Invalid Inputs');
                
            end
        end

        % Button pushed function: SAVEASButton
        function SAVEASButtonPushed(app, event)
            create_outputfile(app);
        end

        % Button pushed function: CLOSEButton
        function CLOSEButtonPushed(app, event)
            delete(app)
        end

        % Button pushed function: RESETButton
        function RESETButtonPushed(app, event)
            
            app.DurationofanimationsecsEditField.Value = 0;
            app.Initialanglemadebycrank1withhorizontaldegreesEditField.Value =0;       
            app.AngularvelocityofDrivercrankradsEditField.Value=0;
            app.InclinationoffixedlinkdegreesEditField.Value = 0;
            app.TimestepsecsEditField.Value=0;
            app.LENGTHmEditField_4.Value=0;
            app.LENGTHmEditField_3.Value=0;
            app.LENGTHmEditField_2.Value=0;
            app.LENGTHmEditField.Value=0;
            app.LENGTHmEditField_5.Value=0;
            app.LENGTHmEditField_6.Value=0;            
         app.CRANKROCKERCheckBox.Value = false;
         app.DRAGLINKCheckBox.Value=false;
            app.NetaccelerationDropDown.Value=app.NetaccelerationDropDown.Items(1);
            app.NetvelocityDropDown.Value=app.NetvelocityDropDown.Items(1);            
            app.AnimationCheckBox.Value = false;
            app.VelocityandaccelerationanalysisCheckBox.Value = false;
            app.StaticforceanalysisCheckBox.Value = false;
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 890 521];
            app.UIFigure.Name = 'MATLAB App';

            % Create INPUTLabel
            app.INPUTLabel = uilabel(app.UIFigure);
            app.INPUTLabel.Position = [107 489 140 32];
            app.INPUTLabel.Text = 'INPUT';

            % Create aLabel
            app.aLabel = uilabel(app.UIFigure);
            app.aLabel.Position = [30 451 101 30];
            app.aLabel.Text = 'a';

            % Create LENGTHmEditFieldLabel
            app.LENGTHmEditFieldLabel = uilabel(app.UIFigure);
            app.LENGTHmEditFieldLabel.HorizontalAlignment = 'right';
            app.LENGTHmEditFieldLabel.Position = [53 455 71 22];
            app.LENGTHmEditFieldLabel.Text = 'LENGTH(m)';

            % Create LENGTHmEditField
            app.LENGTHmEditField = uieditfield(app.UIFigure, 'numeric');
            app.LENGTHmEditField.Position = [133 451 69 29];

            % Create bLabel
            app.bLabel = uilabel(app.UIFigure);
            app.bLabel.Position = [30 401 58 23];
            app.bLabel.Text = 'b';

            % Create LENGTHmEditField_2Label
            app.LENGTHmEditField_2Label = uilabel(app.UIFigure);
            app.LENGTHmEditField_2Label.HorizontalAlignment = 'right';
            app.LENGTHmEditField_2Label.Position = [53 402 71 22];
            app.LENGTHmEditField_2Label.Text = 'LENGTH(m)';

            % Create LENGTHmEditField_2
            app.LENGTHmEditField_2 = uieditfield(app.UIFigure, 'numeric');
            app.LENGTHmEditField_2.Position = [133 398 69 29];

            % Create cLabel
            app.cLabel = uilabel(app.UIFigure);
            app.cLabel.Position = [30 353 58 23];
            app.cLabel.Text = 'c';

            % Create LENGTHmEditField_3Label
            app.LENGTHmEditField_3Label = uilabel(app.UIFigure);
            app.LENGTHmEditField_3Label.HorizontalAlignment = 'right';
            app.LENGTHmEditField_3Label.Position = [53 350 71 22];
            app.LENGTHmEditField_3Label.Text = 'LENGTH(m)';

            % Create LENGTHmEditField_3
            app.LENGTHmEditField_3 = uieditfield(app.UIFigure, 'numeric');
            app.LENGTHmEditField_3.Position = [133 347 69 29];

            % Create dLabel
            app.dLabel = uilabel(app.UIFigure);
            app.dLabel.Position = [30 305 58 23];
            app.dLabel.Text = 'd';

            % Create LENGTHmEditField_4Label
            app.LENGTHmEditField_4Label = uilabel(app.UIFigure);
            app.LENGTHmEditField_4Label.HorizontalAlignment = 'right';
            app.LENGTHmEditField_4Label.Position = [50 302 71 22];
            app.LENGTHmEditField_4Label.Text = 'LENGTH(m)';

            % Create LENGTHmEditField_4
            app.LENGTHmEditField_4 = uieditfield(app.UIFigure, 'numeric');
            app.LENGTHmEditField_4.Position = [130 299 69 29];

            % Create crLabel
            app.crLabel = uilabel(app.UIFigure);
            app.crLabel.Position = [30 254 58 23];
            app.crLabel.Text = 'cr';

            % Create LENGTHmEditField_5Label
            app.LENGTHmEditField_5Label = uilabel(app.UIFigure);
            app.LENGTHmEditField_5Label.HorizontalAlignment = 'right';
            app.LENGTHmEditField_5Label.Position = [53 251 71 22];
            app.LENGTHmEditField_5Label.Text = 'LENGTH(m)';

            % Create LENGTHmEditField_5
            app.LENGTHmEditField_5 = uieditfield(app.UIFigure, 'numeric');
            app.LENGTHmEditField_5.Position = [133 248 69 29];

            % Create ELabel
            app.ELabel = uilabel(app.UIFigure);
            app.ELabel.Position = [30 199 58 23];
            app.ELabel.Text = 'E';

            % Create LENGTHmEditField_6Label
            app.LENGTHmEditField_6Label = uilabel(app.UIFigure);
            app.LENGTHmEditField_6Label.HorizontalAlignment = 'right';
            app.LENGTHmEditField_6Label.Position = [51 202 71 22];
            app.LENGTHmEditField_6Label.Text = 'LENGTH(m)';

            % Create LENGTHmEditField_6
            app.LENGTHmEditField_6 = uieditfield(app.UIFigure, 'numeric');
            app.LENGTHmEditField_6.Position = [131 199 69 29];

            % Create TIMEINTERVALLabel
            app.TIMEINTERVALLabel = uilabel(app.UIFigure);
            app.TIMEINTERVALLabel.Position = [340 480 96 40];
            app.TIMEINTERVALLabel.Text = 'TIME INTERVAL';

            % Create ANGLEVARIATIONLabel
            app.ANGLEVARIATIONLabel = uilabel(app.UIFigure);
            app.ANGLEVARIATIONLabel.Position = [367 336 112 40];
            app.ANGLEVARIATIONLabel.Text = 'ANGLE VARIATION';

            % Create VelocityandaccelerationanalysisCheckBox
            app.VelocityandaccelerationanalysisCheckBox = uicheckbox(app.UIFigure);
            app.VelocityandaccelerationanalysisCheckBox.Text = 'Velocity and acceleration analysis';
            app.VelocityandaccelerationanalysisCheckBox.Position = [666 442 203 32];

            % Create StaticforceanalysisCheckBox
            app.StaticforceanalysisCheckBox = uicheckbox(app.UIFigure);
            app.StaticforceanalysisCheckBox.Text = 'Static force analysis';
            app.StaticforceanalysisCheckBox.Position = [666 411 188 32];

            % Create CALCULATEButton
            app.CALCULATEButton = uibutton(app.UIFigure, 'push');
            app.CALCULATEButton.ButtonPushedFcn = createCallbackFcn(app, @CALCULATEButtonPushed, true);
            app.CALCULATEButton.Position = [682 264 137 31];
            app.CALCULATEButton.Text = 'CALCULATE';

            % Create SAVEASButton
            app.SAVEASButton = uibutton(app.UIFigure, 'push');
            app.SAVEASButton.ButtonPushedFcn = createCallbackFcn(app, @SAVEASButtonPushed, true);
            app.SAVEASButton.Position = [679 160 143 35];
            app.SAVEASButton.Text = 'SAVE AS';

            % Create RESETButton
            app.RESETButton = uibutton(app.UIFigure, 'push');
            app.RESETButton.ButtonPushedFcn = createCallbackFcn(app, @RESETButtonPushed, true);
            app.RESETButton.Position = [679 43 149 35];
            app.RESETButton.Text = 'RESET';

            % Create CLOSEButton
            app.CLOSEButton = uibutton(app.UIFigure, 'push');
            app.CLOSEButton.ButtonPushedFcn = createCallbackFcn(app, @CLOSEButtonPushed, true);
            app.CLOSEButton.Position = [677 101 142 35];
            app.CLOSEButton.Text = 'CLOSE';

            % Create AnimationCheckBox
            app.AnimationCheckBox = uicheckbox(app.UIFigure);
            app.AnimationCheckBox.Text = 'Animation';
            app.AnimationCheckBox.Position = [666 380 188 32];

            % Create TimestepsecsEditFieldLabel
            app.TimestepsecsEditFieldLabel = uilabel(app.UIFigure);
            app.TimestepsecsEditFieldLabel.HorizontalAlignment = 'right';
            app.TimestepsecsEditFieldLabel.Position = [329 398 91 22];
            app.TimestepsecsEditFieldLabel.Text = 'Time step(secs)';

            % Create TimestepsecsEditField
            app.TimestepsecsEditField = uieditfield(app.UIFigure, 'numeric');
            app.TimestepsecsEditField.Position = [435 394 48 30];

            % Create DurationofanimationsecsEditFieldLabel
            app.DurationofanimationsecsEditFieldLabel = uilabel(app.UIFigure);
            app.DurationofanimationsecsEditFieldLabel.HorizontalAlignment = 'right';
            app.DurationofanimationsecsEditFieldLabel.Position = [268 435 152 22];
            app.DurationofanimationsecsEditFieldLabel.Text = 'Duration of animation(secs)';

            % Create DurationofanimationsecsEditField
            app.DurationofanimationsecsEditField = uieditfield(app.UIFigure, 'numeric');
            app.DurationofanimationsecsEditField.Position = [435 431 48 30];

            % Create Initialanglemadebycrank1withhorizontaldegreesEditFieldLabel
            app.Initialanglemadebycrank1withhorizontaldegreesEditFieldLabel = uilabel(app.UIFigure);
            app.Initialanglemadebycrank1withhorizontaldegreesEditFieldLabel.HorizontalAlignment = 'right';
            app.Initialanglemadebycrank1withhorizontaldegreesEditFieldLabel.Position = [219 305 290 22];
            app.Initialanglemadebycrank1withhorizontaldegreesEditFieldLabel.Text = 'Initial angle made by crank1 with horizontal(degrees)';

            % Create Initialanglemadebycrank1withhorizontaldegreesEditField
            app.Initialanglemadebycrank1withhorizontaldegreesEditField = uieditfield(app.UIFigure, 'numeric');
            app.Initialanglemadebycrank1withhorizontaldegreesEditField.Position = [524 301 48 30];

            % Create InclinationoffixedlinkdegreesEditFieldLabel
            app.InclinationoffixedlinkdegreesEditFieldLabel = uilabel(app.UIFigure);
            app.InclinationoffixedlinkdegreesEditFieldLabel.HorizontalAlignment = 'right';
            app.InclinationoffixedlinkdegreesEditFieldLabel.Position = [288 243 174 22];
            app.InclinationoffixedlinkdegreesEditFieldLabel.Text = 'Inclination of fixed link(degrees)';

            % Create InclinationoffixedlinkdegreesEditField
            app.InclinationoffixedlinkdegreesEditField = uieditfield(app.UIFigure, 'numeric');
            app.InclinationoffixedlinkdegreesEditField.Position = [477 239 48 30];

            % Create AngularvelocityofDrivercrankradsEditFieldLabel_2
            app.AngularvelocityofDrivercrankradsEditFieldLabel_2 = uilabel(app.UIFigure);
            app.AngularvelocityofDrivercrankradsEditFieldLabel_2.HorizontalAlignment = 'right';
            app.AngularvelocityofDrivercrankradsEditFieldLabel_2.Position = [277 202 204 29];
            app.AngularvelocityofDrivercrankradsEditFieldLabel_2.Text = 'Angular velocity of Driver crank(rad/s)';

            % Create AngularvelocityofDrivercrankradsEditField
            app.AngularvelocityofDrivercrankradsEditField = uieditfield(app.UIFigure, 'numeric');
            app.AngularvelocityofDrivercrankradsEditField.Position = [480 203 45 30];

            % Create ANALYSISLabel
            app.ANALYSISLabel = uilabel(app.UIFigure);
            app.ANALYSISLabel.Position = [718 473 60 37];
            app.ANALYSISLabel.Text = 'ANALYSIS';

            % Create NetvelocityDropDownLabel
            app.NetvelocityDropDownLabel = uilabel(app.UIFigure);
            app.NetvelocityDropDownLabel.HorizontalAlignment = 'right';
            app.NetvelocityDropDownLabel.Position = [666 345 67 22];
            app.NetvelocityDropDownLabel.Text = 'Net velocity';

            % Create NetvelocityDropDown
            app.NetvelocityDropDown = uidropdown(app.UIFigure);
            app.NetvelocityDropDown.Items = {'Select', 'Point A', 'Point B', 'Slider'};
            app.NetvelocityDropDown.Position = [748 345 106 22];
            app.NetvelocityDropDown.Value = 'Select';

            % Create NetaccelerationDropDownLabel
            app.NetaccelerationDropDownLabel = uilabel(app.UIFigure);
            app.NetaccelerationDropDownLabel.HorizontalAlignment = 'right';
            app.NetaccelerationDropDownLabel.Position = [644 315 92 22];
            app.NetaccelerationDropDownLabel.Text = 'Net acceleration';

            % Create NetaccelerationDropDown
            app.NetaccelerationDropDown = uidropdown(app.UIFigure);
            app.NetaccelerationDropDown.Items = {'Select', 'Point A', 'PointB', 'Slider'};
            app.NetaccelerationDropDown.Position = [751 315 106 22];
            app.NetaccelerationDropDown.Value = 'Select';

            % Create asEditFieldLabel
            app.asEditFieldLabel = uilabel(app.UIFigure);
            app.asEditFieldLabel.HorizontalAlignment = 'right';
            app.asEditFieldLabel.Position = [258 150 25 22];
            app.asEditFieldLabel.Text = 'as';

            % Create asEditField
            app.asEditField = uieditfield(app.UIFigure, 'numeric');
            app.asEditField.Position = [288 147 42 27];

            % Create bsEditFieldLabel
            app.bsEditFieldLabel = uilabel(app.UIFigure);
            app.bsEditFieldLabel.HorizontalAlignment = 'right';
            app.bsEditFieldLabel.Position = [258 113 25 22];
            app.bsEditFieldLabel.Text = 'bs';

            % Create bsEditField
            app.bsEditField = uieditfield(app.UIFigure, 'numeric');
            app.bsEditField.Position = [287 110 38 27];

            % Create csEditFieldLabel
            app.csEditFieldLabel = uilabel(app.UIFigure);
            app.csEditFieldLabel.HorizontalAlignment = 'right';
            app.csEditFieldLabel.Position = [258 80 25 22];
            app.csEditFieldLabel.Text = 'cs';

            % Create csEditField
            app.csEditField = uieditfield(app.UIFigure, 'numeric');
            app.csEditField.Position = [287 77 40 27];

            % Create crsEditFieldLabel
            app.crsEditFieldLabel = uilabel(app.UIFigure);
            app.crsEditFieldLabel.HorizontalAlignment = 'right';
            app.crsEditFieldLabel.Position = [246 46 33 22];
            app.crsEditFieldLabel.Text = 'crs';

            % Create crsEditField
            app.crsEditField = uieditfield(app.UIFigure, 'numeric');
            app.crsEditField.Position = [284 43 42 27];

            % Create FasEditFieldLabel
            app.FasEditFieldLabel = uilabel(app.UIFigure);
            app.FasEditFieldLabel.HorizontalAlignment = 'right';
            app.FasEditFieldLabel.Position = [340 161 26 22];
            app.FasEditFieldLabel.Text = 'Fas';

            % Create FasEditField
            app.FasEditField = uieditfield(app.UIFigure, 'numeric');
            app.FasEditField.Position = [368 158 50 27];

            % Create FbsEditFieldLabel
            app.FbsEditFieldLabel = uilabel(app.UIFigure);
            app.FbsEditFieldLabel.HorizontalAlignment = 'right';
            app.FbsEditFieldLabel.Position = [339 124 26 22];
            app.FbsEditFieldLabel.Text = 'Fbs';

            % Create FbsEditField
            app.FbsEditField = uieditfield(app.UIFigure, 'numeric');
            app.FbsEditField.Position = [367 121 50 27];

            % Create FcsEditFieldLabel
            app.FcsEditFieldLabel = uilabel(app.UIFigure);
            app.FcsEditFieldLabel.HorizontalAlignment = 'right';
            app.FcsEditFieldLabel.Position = [338 87 25 22];
            app.FcsEditFieldLabel.Text = 'Fcs';

            % Create FcsEditField
            app.FcsEditField = uieditfield(app.UIFigure, 'numeric');
            app.FcsEditField.Position = [365 84 49 27];

            % Create FcrsEditFieldLabel
            app.FcrsEditFieldLabel = uilabel(app.UIFigure);
            app.FcrsEditFieldLabel.HorizontalAlignment = 'right';
            app.FcrsEditFieldLabel.Position = [334 55 29 22];
            app.FcrsEditFieldLabel.Text = 'Fcrs';

            % Create FcrsEditField
            app.FcrsEditField = uieditfield(app.UIFigure, 'numeric');
            app.FcrsEditField.Position = [366 52 49 27];

            % Create FssEditFieldLabel
            app.FssEditFieldLabel = uilabel(app.UIFigure);
            app.FssEditFieldLabel.HorizontalAlignment = 'right';
            app.FssEditFieldLabel.Position = [339 21 25 22];
            app.FssEditFieldLabel.Text = 'Fss';

            % Create FssEditField
            app.FssEditField = uieditfield(app.UIFigure, 'numeric');
            app.FssEditField.Position = [367 18 49 27];

            % Create thetaFasEditFieldLabel
            app.thetaFasEditFieldLabel = uilabel(app.UIFigure);
            app.thetaFasEditFieldLabel.HorizontalAlignment = 'right';
            app.thetaFasEditFieldLabel.Position = [427 161 52 22];
            app.thetaFasEditFieldLabel.Text = 'thetaFas';

            % Create thetaFasEditField
            app.thetaFasEditField = uieditfield(app.UIFigure, 'numeric');
            app.thetaFasEditField.Position = [478 158 47 27];

            % Create thetaFbsEditFieldLabel
            app.thetaFbsEditFieldLabel = uilabel(app.UIFigure);
            app.thetaFbsEditFieldLabel.HorizontalAlignment = 'right';
            app.thetaFbsEditFieldLabel.Position = [427 129 52 22];
            app.thetaFbsEditFieldLabel.Text = 'thetaFbs';

            % Create thetaFbsEditField
            app.thetaFbsEditField = uieditfield(app.UIFigure, 'numeric');
            app.thetaFbsEditField.Position = [480 126 45 27];

            % Create thetaFcsEditFieldLabel
            app.thetaFcsEditFieldLabel = uilabel(app.UIFigure);
            app.thetaFcsEditFieldLabel.HorizontalAlignment = 'right';
            app.thetaFcsEditFieldLabel.Position = [429 98 51 22];
            app.thetaFcsEditFieldLabel.Text = 'thetaFcs';

            % Create thetaFcsEditField
            app.thetaFcsEditField = uieditfield(app.UIFigure, 'numeric');
            app.thetaFcsEditField.Position = [481 95 44 27];

            % Create thetaFcrsEditFieldLabel
            app.thetaFcrsEditFieldLabel = uilabel(app.UIFigure);
            app.thetaFcrsEditFieldLabel.HorizontalAlignment = 'right';
            app.thetaFcrsEditFieldLabel.Position = [427 61 55 22];
            app.thetaFcrsEditFieldLabel.Text = 'thetaFcrs';

            % Create thetaFcrsEditField
            app.thetaFcrsEditField = uieditfield(app.UIFigure, 'numeric');
            app.thetaFcrsEditField.Position = [481 58 46 27];

            % Create thetaFssEditFieldLabel
            app.thetaFssEditFieldLabel = uilabel(app.UIFigure);
            app.thetaFssEditFieldLabel.HorizontalAlignment = 'right';
            app.thetaFssEditFieldLabel.Position = [427 29 51 22];
            app.thetaFssEditFieldLabel.Text = 'thetaFss';

            % Create thetaFssEditField
            app.thetaFssEditField = uieditfield(app.UIFigure, 'numeric');
            app.thetaFssEditField.Position = [479 26 44 27];

            % Create TYPEFORLOOP1Label
            app.TYPEFORLOOP1Label = uilabel(app.UIFigure);
            app.TYPEFORLOOP1Label.Position = [50 123 115 22];
            app.TYPEFORLOOP1Label.Text = 'TYPE FOR LOOP1';

            % Create CRANKROCKERCheckBox
            app.CRANKROCKERCheckBox = uicheckbox(app.UIFigure);
            app.CRANKROCKERCheckBox.Text = 'CRANK ROCKER';
            app.CRANKROCKERCheckBox.Position = [54 92 119 22];

            % Create DRAGLINKCheckBox
            app.DRAGLINKCheckBox = uicheckbox(app.UIFigure);
            app.DRAGLINKCheckBox.Text = 'DRAG LINK';
            app.DRAGLINKCheckBox.Position = [54 67 87 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Variable_stroke_mechanism_t2_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
classdef Variable_Stroke_mechanism_t1_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        thetaFc1sEditField            matlab.ui.control.NumericEditField
        thetaFc1sEditFieldLabel       matlab.ui.control.Label
        thetaFb1sEditField            matlab.ui.control.NumericEditField
        thetaFb1sEditFieldLabel       matlab.ui.control.Label
        thetaFcsEditField             matlab.ui.control.NumericEditField
        thetaFcsEditFieldLabel        matlab.ui.control.Label
        thetaFbsEditField             matlab.ui.control.NumericEditField
        thetaFbsEditFieldLabel        matlab.ui.control.Label
        thetaFasEditField             matlab.ui.control.NumericEditField
        thetaFasEditFieldLabel        matlab.ui.control.Label
        c1sEditField                  matlab.ui.control.NumericEditField
        c1sEditFieldLabel             matlab.ui.control.Label
        b1sEditField                  matlab.ui.control.NumericEditField
        b1sEditFieldLabel             matlab.ui.control.Label
        csEditField                   matlab.ui.control.NumericEditField
        csEditFieldLabel              matlab.ui.control.Label
        bsEditField                   matlab.ui.control.NumericEditField
        bsEditFieldLabel              matlab.ui.control.Label
        asEditField                   matlab.ui.control.NumericEditField
        asEditFieldLabel              matlab.ui.control.Label
        Fc1sEditField                 matlab.ui.control.NumericEditField
        Fc1sEditFieldLabel            matlab.ui.control.Label
        Fb1sEditField                 matlab.ui.control.NumericEditField
        Fb1sEditFieldLabel            matlab.ui.control.Label
        FcsEditField                  matlab.ui.control.NumericEditField
        FcsEditFieldLabel             matlab.ui.control.Label
        FbsEditField                  matlab.ui.control.NumericEditField
        FbsEditFieldLabel             matlab.ui.control.Label
        FasEditField                  matlab.ui.control.NumericEditField
        FasEditFieldLabel             matlab.ui.control.Label
        NetaccelerationDropDown       matlab.ui.control.DropDown
        NetaccelerationDropDownLabel  matlab.ui.control.Label
        NetvelocityDropDown           matlab.ui.control.DropDown
        NetvelocityDropDownLabel      matlab.ui.control.Label
        AngularaccelarationDropDown   matlab.ui.control.DropDown
        AngularaccelarationDropDownLabel  matlab.ui.control.Label
        AngularvelocityDropDown       matlab.ui.control.DropDown
        AngularvelocityDropDownLabel  matlab.ui.control.Label
        PLOTSLabel                    matlab.ui.control.Label
        Initialanglemadebycrank1withhorizontaldegreesEditField  matlab.ui.control.NumericEditField
        Initialanglemadebycrank1withhorizontaldegreesEditFieldLabel  matlab.ui.control.Label
        AnimationCheckBox             matlab.ui.control.CheckBox
        Inclinationoffixedlink2degreesEditField  matlab.ui.control.NumericEditField
        Inclinationoffixedlink2degreesEditFieldLabel  matlab.ui.control.Label
        LENGTHmEditField_7            matlab.ui.control.NumericEditField
        LENGTHmEditField_7Label       matlab.ui.control.Label
        LENGTHmEditField_6            matlab.ui.control.NumericEditField
        LENGTHmEditField_6Label       matlab.ui.control.Label
        LENGTHmEditField_5            matlab.ui.control.NumericEditField
        LENGTHmEditField_5Label       matlab.ui.control.Label
        d1Label                       matlab.ui.control.Label
        c1Label                       matlab.ui.control.Label
        b1Label                       matlab.ui.control.Label
        AngularvelocityofDrivercrankradsEditField  matlab.ui.control.NumericEditField
        AngularvelocityofDrivercrankradsEditFieldLabel  matlab.ui.control.Label
        CLOSEButton                   matlab.ui.control.Button
        RESETButton                   matlab.ui.control.Button
        SAVEASButton                  matlab.ui.control.Button
        CALCULATEButton               matlab.ui.control.Button
        StaticforceanalysisCheckBox   matlab.ui.control.CheckBox
        VelocityandaccelerationanalysisCheckBox  matlab.ui.control.CheckBox
        Inclinationoffixedlink1degreesEditField  matlab.ui.control.NumericEditField
        Inclinationoffixedlink1degreesEditFieldLabel  matlab.ui.control.Label
        ANGLEVARIATIONLabel           matlab.ui.control.Label
        TimestepsecsEditField         matlab.ui.control.NumericEditField
        TimestepsecsEditFieldLabel    matlab.ui.control.Label
        DurationofanimationsecsEditField  matlab.ui.control.NumericEditField
        DurationofanimationsecsEditFieldLabel  matlab.ui.control.Label
        TIMEINTERVALLabel             matlab.ui.control.Label
        ANALYSISLabel                 matlab.ui.control.Label
        LENGTHmEditField_4            matlab.ui.control.NumericEditField
        LENGTHmEditField_4Label       matlab.ui.control.Label
        dLabel                        matlab.ui.control.Label
        LENGTHmEditField_3            matlab.ui.control.NumericEditField
        LENGTHmEditField_3Label       matlab.ui.control.Label
        cLabel                        matlab.ui.control.Label
        LENGTHmEditField_2            matlab.ui.control.NumericEditField
        LENGTHmEditField_2Label       matlab.ui.control.Label
        bLabel                        matlab.ui.control.Label
        aLabel                        matlab.ui.control.Label
        LENGTHmEditField_8            matlab.ui.control.NumericEditField
        LENGTHmEditField_8Label       matlab.ui.control.Label
        INPUTLabel                    matlab.ui.control.Label
    end

methods (Access = public)
    function Animation_t1(app)
            l1 = app.LENGTHmEditField_8.Value;
            l2 = app.LENGTHmEditField_2.Value;
            l3= app.LENGTHmEditField_3.Value;
            l4= app.LENGTHmEditField_4.Value;
            
            l5 = app.LENGTHmEditField_5.Value;
            l6 = app.LENGTHmEditField_6.Value;
            l7 = app.LENGTHmEditField_7.Value;
            
            B=[l3;l5;l6;l7];
            Q=sort(B);
            
            A = [l1;l2;l3;l4];            
            P=sort(A);
            
            if(P(1)+P(4)<P(2)+P(3))
               a=P(1);
               b=P(4);
               c=P(2);
               d=P(3);
               if(Q(1)+Q(4)<Q(2)+Q(3))
                 if(c==Q(4))
                   a1=Q(4);
                   b1=Q(2);
                   c1=Q(3);
                   d1=Q(1); 
                 else 
                    b1=Q(1);        
                    c1=Q(4);
                    if(P(2)==Q(2))
                        a1=Q(2);
                        d1 = Q(3);
                    else d1 = Q(2);
                         a1 = Q(3);
                    end 
                 end
                 check_3=1;
               elseif(Q(1)+Q(4)>Q(2)+Q(3))
                   if(c==Q(1)) 
                        b1=Q(4)
                        c1=Q(2)
                        d1=Q(2)
                        a1=Q(1)
                        check_3=1;
                    else 
                      d1=Q(1);
                      a1=P(2);
                      b1=Q(4);
                      c1=Q(3);
                      check_3=1;
                    end                                                         
               else check_3 = 0;
                    
               end
            else check_3 =0;
                                 
            end
        
            
            t = app.DurationofanimationsecsEditField.Value;
            delt =app.TimestepsecsEditField.Value;
            
            initial_angle = app.Initialanglemadebycrank1withhorizontaldegreesEditField.Value;
            Wa = app.AngularvelocityofDrivercrankradsEditField.Value;
            theta1 = app.Inclinationoffixedlink1degreesEditField.Value;
            theta2 = app.Inclinationoffixedlink2degreesEditField.Value;
            
            A_a=0;
            
            final_angle = (((Wa*t*180/pi)+initial_angle));  
            
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
            
            k1=((a1^2-b1^2+c1^2+d1^2)/2);
            A1=-a1.*(d1.*cosd(theta2)-c1).*cosd(si)+k1-(c1*d1*cosd(theta2))-a1*d1.*sind(si).*sind(theta2);
            B1=-(2.*a1.*c1.*sind(si)-2*c1*d1.*sind(theta2));
            C1=-a1*(d1*cosd(theta2)+c1).*cosd(si)+k1+(c1*d1.*cosd(theta2))-a1.*d1.*sind(si).*sind(theta2);
            si1=2.*atand((-B1-sqrt((B1.^2)-4.*A1.*C1))./(2*A1));
            
            G1=((a1^2+b1^2-c1^2+d1^2)/2);
            D1=-a1.*(d1.*cosd(theta2)+b1).*cosd(si)+G1+(d1*cosd(theta2)*b1)-a1*d1*sind(theta2).*sind(si);
            E1=2*a1*b1.*sind(si)-2*b1*d1.*sind(theta2);
            F1=-(a1.*(d1*cosd(theta2)-b1).*cosd(si))+G1-(b1*d1*cosd(theta2))-a1.*d1.*sind(si).*sind(theta2);
            beta1=2.*atand((-E1+sqrt((E1.^2)-4.*D1.*F1))./(2*D1));
            
            
            
             plot([0 a*cosd(theta)], [0 a*sind(theta)],'ko-','LineWidth',2);hold on;
             plot([d*cosd(theta1) d*cosd(theta1)+c*cosd(si)], [d*sind(theta1) d*sind(theta1)+c*sind(si)], 'bo-','LineWidth',2); hold on;
             plot([0 d*cosd(theta1)], [0 d*sind(theta1)], 'mo-','LineWidth',2); hold on;
             plot([a*cosd(theta) a*cosd(theta)+b*cosd(beta)], [a*sind(theta) a*sind(theta)+b*sind(beta)], 'ro-','LineWidth',2);hold on;
             plot([d*cosd(theta1) d*cosd(theta1)+d1*cosd(theta2)], [d*sind(theta1) d*sind(theta1)+d1*sind(theta2)], 'go-','LineWidth',2);hold on;
             plot([d*cosd(theta1)+d1*cosd(theta2) d*cosd(theta1)+d1*cosd(theta2)+c1*cosd(si1)],[d*sind(theta1)+d1*sind(theta2) d*sind(theta1)+d1*sind(theta2)+c1*sind(si1)], 'co-','LineWidth',2);hold on;
             plot([d*cosd(theta1)+c*cosd(si) d*cosd(theta1)+c*cosd(si)+b1*cosd(beta1)], [d*sind(theta1)+c*sind(si) d*sind(theta1)+c*sind(si)+b1*sind(beta1)], 'yo-','LineWidth',2);hold off;
             grid on
             axis([-2*a d1+c1+d+a -2*a d1+c1+d+a]);
             pbaspect([1 1 1]);
             pause(0.01);
             drawnow
         end
            
        end
    
        
        function Velocity_accleration_A(app)
                        
            l1 = app.LENGTHmEditField_8.Value;
            l2 = app.LENGTHmEditField_2.Value;
            l3= app.LENGTHmEditField_3.Value;
            l4= app.LENGTHmEditField_4.Value;
            
            l5 = app.LENGTHmEditField_5.Value;
            l6 = app.LENGTHmEditField_6.Value;
            l7 = app.LENGTHmEditField_7.Value;
            
            B=[l3;l5;l6;l7];
            Q=sort(B);
            
            A = [l1;l2;l3;l4];            
            P=sort(A);
            
            if(P(1)+P(4)<P(2)+P(3))
               a=P(1);
               b=P(4);
               c=P(2);
               d=P(3);
               if(Q(1)+Q(4)<Q(2)+Q(3))
                 if(c==Q(4))
                   a1=Q(4);
                   b1=Q(2);
                   c1=Q(3);
                   d1=Q(1); 
                 else 
                    b1=Q(1);        
                    c1=Q(4);
                    if(P(2)==Q(2))
                        a1=Q(2);
                        d1 = Q(3);
                    else d1 = Q(2);
                         a1 = Q(3);
                    end 
                 end
                 check_3=1;
               elseif(Q(1)+Q(4)>Q(2)+Q(3))
                 if(c==Q(1)) 
                        b1=Q(4)
                        c1=Q(2)
                        d1=Q(2)
                        a1=Q(1)
                        check_3=1;
                    else 
                      d1=Q(1);
                      a1=P(2);
                      b1=Q(4);
                      c1=Q(3);
                      check_3=1;
                 end                   
               else check_3 = 0;
                    
               end
            else check_3 =0;
                                 
            end
            
            time_of_running = app.DurationofanimationsecsEditField.Value;
            delt =app.TimestepsecsEditField.Value;
            
            Wa = app.AngularvelocityofDrivercrankradsEditField.Value;
            theta1 = app.Inclinationoffixedlink1degreesEditField.Value;
            theta2 = app.Inclinationoffixedlink2degreesEditField.Value;
            initial_angle = app.Initialanglemadebycrank1withhorizontaldegreesEditField.Value;

            Aa=0;           
            
            final_angle=initial_angle+((time_of_running*Wa)*180/pi);
            
            theta=initial_angle:(Wa*180/pi)/2:final_angle;            
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
            
            a1=c;
            k1=((a1^2-b1^2+c1^2+d1^2)/2);
            A1=-a1.*(d1.*cosd(theta2)-c1).*cosd(si)+k1-(c1*d1*cosd(theta2))-a1*d1.*sind(si).*sind(theta2);
            B1=-(2.*a1.*c1.*sind(si)-2*c1*d1.*sind(theta2));
            C1=-a1*(d1*cosd(theta2)+c1).*cosd(si)+k1+(c1*d1.*cosd(theta2))-a1.*d1.*sind(si).*sind(theta2);
            si1=2.*atand((-B1-sqrt((B1.^2)-4.*A1.*C1))./(2*A1));
            
            G1=((a1^2+b1^2-c1^2+d1^2)/2);
            D1=-a1.*(d1.*cosd(theta2)+b1).*cosd(si)+G1+(d1*cosd(theta2)*b1)-a1*d1*sind(theta2).*sind(si);
            E1=2*a1*b1.*sind(si)-2*b1*d1.*sind(theta2);
            F1=-(a1.*(d1*cosd(theta2)-b1).*cosd(si))+G1-(b1*d1*cosd(theta2))-a1.*d1.*sind(si).*sind(theta2);
            beta1=2.*atand((-E1+sqrt((E1.^2)-4.*D1.*F1))./(2*D1));
            
            Wc=(a.*Wa.*sind(beta-theta))./(c.*sind(beta-si));
            Wb=(-a*Wa.*sind(si-theta))./(b.*sind(si-beta));
            Ab=((a.*Aa.*sind(si-theta))-(a.*(Wa.^2).*cosd(si-theta))-(b.*(Wb.^2).*cosd(si-beta))-(c.*(Wc.^2)))./(b.*sind(beta-si));
            Ac=((a.*Aa.*sind(beta-theta))-(a.*(Wa.^2).*cosd(beta-theta))-(b.*(Wb.^2))+(c.*(Wc.^2).*cosd(beta-si)))./(c.*sind(beta-si));
            Wc1=(a1.*Wc.*sind(beta1-si))./(c1.*sind(beta1-si1));
            Wb1=(-a1*Wc.*sind(si1-si))./(b1.*sind(si1-beta1));
            Ab1=((a1.*Ac.*sind(si1-si))-(a1.*(Wc.^2).*cosd(si1-si))-(b1.*(Wb1.^2).*cosd(si1-beta1))-(c1.*(Wc1.^2)))./(b1.*sind(beta1-si1));
            Ac1=((a1.*Ac.*sind(beta1-si))-(a1.*(Wc.^2).*cosd(beta1-si))-(b1.*(Wb1.^2))+(c1.*(Wc1.^2).*cosd(beta1-si1)))./(c1.*sind(beta1-si1));
            

            
            %positions of each point
            O=[0;0];
            A=[a*cosd(theta);a*sind(theta)];
            B=[d*cosd(theta1)+c*cosd(si);d*sind(theta1)+c*sind(si)];
            C=[d*cosd(theta1);d*sind(theta1)];
            D=[d*cosd(theta1)+c*cosd(si)+b1*cosd(beta1);d*sind(theta1)+c*sind(si)+b1*sind(beta1)];
            E=[d*cosd(theta1)+d1*cosd(theta2);d*sind(theta1)+d1*sind(theta2)];
            
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
                      
                       
            D_x=D(1,:);
            D_y=D(2,:);
            
            D_vx=diff(D_x)./diff(time);
            D_vy=diff(D_y)./diff(time);
            
            D_v=sqrt(D_vx.^2+D_vy.^2);
            D_v=[0 D_v];
            
            D_a=diff(D_v)./diff(time);
            D_a=[0 D_a];            
            
            
           if app.AngularvelocityDropDown.Value
               value = app.AngularvelocityDropDown.Value;
                   if strcmpi(value,'ROCKER1')
                         plot(time,Wc);
                         ylabel('Angular velocity of rocker1 (rad/s)')
                         xlabel('time (secs)')
                         grid on;
                         
                   elseif strcmpi(value,'COUPLER1')
                       plot(time,Wb);grid on
                       ylabel('Angular velocity of coupler1 (rad/s)')
                       xlabel('time (secs)')
                       
                   elseif strcmpi(value,'ROCKER2')                       
                       plot(time,Wc1);grid on
                       ylabel('Angular velocity of rocker2 (rad/s)')
                       xlabel('time (secs)')
                       
                   elseif strcmpi(value,'COUPLER2')                       
                       plot(time,Wb1);grid on
                       ylabel('Angular velocity of coupler2 (rad/s)')
                       xlabel('time (secs)')
                   end   
           end
               
                   
                   
           if app.AngularaccelarationDropDown.Value
               value2 = app.AngularaccelarationDropDown.Value;
                   if strcmpi(value2,'ROCKER1')                          
                         plot(time,Ac);
                         ylabel('Angular acceleration of rocker1 (rad/s^2)')
                         xlabel('time (secs)')
                         grid on;
                   elseif strcmpi(value2,'COUPLER1')                       
                       plot(time,Ab);
                       ylabel('Angular acceleration of coupler1 (rad/s^2)')
                       xlabel('time (secs)')
                       grid on;
                   elseif strcmpi(value2,'ROCKER2')                       
                       plot(time,Ac1);
                       ylabel('Angular acceleration of rocker2 (rad/s^2)')
                       xlabel('time (secs)')
                       grid on;
                   elseif strcmpi(value2,'COUPLER2')                       
                       plot(time,Ab1);
                       ylabel('Angular acceleration of coupler2 (rad/s^2)')
                       xlabel('time (secs)')
                       grid on;
                   end
           end
               
                   
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
                   elseif strcmpi(value3,'Point D')                        
                       plot(time,D_v);grid on
                       ylabel('Velocity of point D (m/s)')
                       xlabel('time (secs)')
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
                   elseif strcmpi(value4,'Point D')                       
                       plot(time,D_a); grid on
                       ylabel('acceleration of point D (m/s^2)')
                       xlabel('time (secs)')
                         
                   end
           end
                          
            
        end
        
        
        function [check] = check_conditions_t1(app)
            
            l1 = app.LENGTHmEditField_8.Value;
            l2 = app.LENGTHmEditField_2.Value;
            l3= app.LENGTHmEditField_3.Value;
            l4= app.LENGTHmEditField_4.Value;
            
            l5 = app.LENGTHmEditField_5.Value;
            l6 = app.LENGTHmEditField_6.Value;
            l7 = app.LENGTHmEditField_7.Value;
            
            B=[l3;l5;l6;l7];
            Q=sort(B);
            
            A = [l1;l2;l3;l4];            
            P=sort(A);
            
            if(P(1)+P(4)<P(2)+P(3))
               a=P(1);
               b=P(4);
               c=P(2);
               d=P(3);
               if(Q(1)+Q(4)<Q(2)+Q(3))
                 if(c==Q(4))
                   a1=Q(4)
                   b1=Q(2)
                   c1=Q(3)
                   d1=Q(1) 
                 else 
                    b1=Q(1)        
                    c1=Q(4)
                    if(P(2)==Q(2))
                        a1=Q(2)
                        d1 = Q(3)
                    else d1 = Q(2)
                         a1 = Q(3)
                    end 
                 end
                 check_1=1;
               elseif(Q(1)+Q(4)>Q(2)+Q(3))    % elseif(q1+q4==q2+q3)
                    if(c==Q(1)) 
                        b1=Q(4)
                        c1=Q(2)
                        d1=Q(2)
                        a1=Q(1)
                        check_1=1;
                    else 
                      d1=Q(1);
                      a1=P(2);
                      b1=Q(4);
                      c1=Q(3);
                      check_1=1;
                    end                                     
                    
               else check_1 = 0;
                    
               end
            else check_1 =0;
                 check = 0;
                 return;
                                              
            end
            
            theta1 = app.Inclinationoffixedlink1degreesEditField.Value;
            theta2 = app.Inclinationoffixedlink2degreesEditField.Value;
            
            Wa = app.AngularvelocityofDrivercrankradsEditField.Value;
            initial_angle = app.Initialanglemadebycrank1withhorizontaldegreesEditField.Value;
            time_of_running = app.DurationofanimationsecsEditField.Value;
            final_angle=initial_angle+((time_of_running*Wa)*180/pi);
            
         if(check_1 ==1)
            
            for theta=initial_angle:Wa/2:final_angle
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
                
                a1 = c;
                k1=((a1^2-b1^2+c1^2+d1^2)/2);
                A1=-a1.*(d1.*cosd(theta2)-c1).*cosd(si)+k1-(c1*d1*cosd(theta2))-a1*d1.*sind(si).*sind(theta2);
                B1=-(2.*a1.*c1.*sind(si)-2*c1*d1.*sind(theta2));
                C1=-a1*(d1*cosd(theta2)+c1).*cosd(si)+k1+(c1*d1.*cosd(theta2))-a1.*d1.*sind(si).*sind(theta2);
                si1=2.*atand((-B1-sqrt((B1.^2)-4.*A1.*C1))./(2*A1));
                
                G1=((a1^2+b1^2-c1^2+d1^2)/2);
                D1=-a1.*(d1.*cosd(theta2)+b1).*cosd(si)+G1+(d1*cosd(theta2)*b1)-a1*d1*sind(theta2).*sind(si);
                E1=2*a1*b1.*sind(si)-2*b1*d1.*sind(theta2);
                F1=-(a1.*(d1*cosd(theta2)-b1).*cosd(si))+G1-(b1*d1*cosd(theta2))-a1.*d1.*sind(si).*sind(theta2);
                beta1=2.*atand((-E1+sqrt((E1.^2)-4.*D1.*F1))./(2*D1));
            end               

         end
            
            if  ((imag(beta) || imag(si) || imag(beta1) || imag(si1)))
                      
                check_2 = (~((imag(beta) || imag(si) || imag(beta1) || imag(si1)) ));
                check = (check_1 && check_2) ;
                           
            else check_2 = 1;  
                check = (check_1 && check_2) ;
            
            end
            
            return;
        end
    
        
        function create_table(app)
            l1 = app.LENGTHmEditField_8.Value;
            l2 = app.LENGTHmEditField_2.Value;
            l3= app.LENGTHmEditField_3.Value;
            l4= app.LENGTHmEditField_4.Value;
            
            l5 = app.LENGTHmEditField_5.Value;
            l6 = app.LENGTHmEditField_6.Value;
            l7 = app.LENGTHmEditField_7.Value;
            
            B=[l3;l5;l6;l7];
            Q=sort(B);
            
            A = [l1;l2;l3;l4];            
            P=sort(A);
            
            if(P(1)+P(4)<P(2)+P(3))
               a=P(1);
               t1 = "CRANK1";
               b=P(4);
               t2 = "COUPLER1";
               c=P(2);
               t3 = "ROCKER1";
               d=P(3);
               t4 = "FIXED LINK1";
               
               if(Q(1)+Q(4)<Q(2)+Q(3))
                 if(c==Q(4))
                   a1=Q(4);
                   b1=Q(2);
                   t5 = "COUPLER2";
                   c1=Q(3);
                   t6 = "ROCKER2";
                   d1=Q(1);
                   t7 = "FIXED LINK2";
                 else 
                    b1=Q(1);
                    t5 = "COUPLER2";
                    c1=Q(4);
                    t6 = "ROCKER2";
                    if(P(2)==Q(2))
                        a1=Q(2);
                        d1 = Q(3);
                        t7 = "FIXED LINK2";
                    else d1 = Q(2);
                         t7 = "FIXED LINK2";
                         a1 = Q(3);
                    end 
                 end
                 check_3=1;
               elseif(Q(1)+Q(4)>Q(2)+Q(3))
                    if(c==Q(1)) 
                        b1=Q(4)
                        t5 = "COUPLER2";
                        c1=Q(2)
                        t6 = "ROCKER2";
                        d1=Q(2)
                        t7 = "FIXED LINK2";
                        a1=Q(1)
                        check_3=1;
                    else 
                      d1=Q(1);
                      t7 = "FIXED LINK2";
                      a1=P(2);
                      b1=Q(4);
                      t5 = "COUPLER2";
                      c1=Q(3);
                      t6 = "ROCKER2";
                      check_3=1;
                    end     
                   
               else check_3 = 0;
                    
               end
            else check_3 =0;
                                 
            end
            type = [t1;t2;t3;t4;t5;t6;t7];
            time_of_running = app.DurationofanimationsecsEditField.Value;
            delt =app.TimestepsecsEditField.Value;
            
            Wa = app.AngularvelocityofDrivercrankradsEditField.Value;
            theta1 = app.Inclinationoffixedlink1degreesEditField.Value;
            theta2 = app.Inclinationoffixedlink2degreesEditField.Value;
            initial_angle = app.Initialanglemadebycrank1withhorizontaldegreesEditField.Value;

            Aa=0;           
            
            final_angle=initial_angle+((time_of_running*Wa)*180/pi);
            
            theta=initial_angle:(Wa*180/pi)/2:final_angle;            
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
            
            a1=c;
            k1=((a1^2-b1^2+c1^2+d1^2)/2);
            A1=-a1.*(d1.*cosd(theta2)-c1).*cosd(si)+k1-(c1*d1*cosd(theta2))-a1*d1.*sind(si).*sind(theta2);
            B1=-(2.*a1.*c1.*sind(si)-2*c1*d1.*sind(theta2));
            C1=-a1*(d1*cosd(theta2)+c1).*cosd(si)+k1+(c1*d1.*cosd(theta2))-a1.*d1.*sind(si).*sind(theta2);
            si1=2.*atand((-B1-sqrt((B1.^2)-4.*A1.*C1))./(2*A1));
            
            G1=((a1^2+b1^2-c1^2+d1^2)/2);
            D1=-a1.*(d1.*cosd(theta2)+b1).*cosd(si)+G1+(d1*cosd(theta2)*b1)-a1*d1*sind(theta2).*sind(si);
            E1=2*a1*b1.*sind(si)-2*b1*d1.*sind(theta2);
            F1=-(a1.*(d1*cosd(theta2)-b1).*cosd(si))+G1-(b1*d1*cosd(theta2))-a1.*d1.*sind(si).*sind(theta2);
            beta1=2.*atand((-E1+sqrt((E1.^2)-4.*D1.*F1))./(2*D1));
            
            Wc=(a.*Wa.*sind(beta-theta))./(c.*sind(beta-si));
            Wb=(-a*Wa.*sind(si-theta))./(b.*sind(si-beta));
            Ab=((a.*Aa.*sind(si-theta))-(a.*(Wa.^2).*cosd(si-theta))-(b.*(Wb.^2).*cosd(si-beta))-(c.*(Wc.^2)))./(b.*sind(beta-si));
            Ac=((a.*Aa.*sind(beta-theta))-(a.*(Wa.^2).*cosd(beta-theta))-(b.*(Wb.^2))+(c.*(Wc.^2).*cosd(beta-si)))./(c.*sind(beta-si));
            Wc1=(a1.*Wc.*sind(beta1-si))./(c1.*sind(beta1-si1));
            Wb1=(-a1*Wc.*sind(si1-si))./(b1.*sind(si1-beta1));
            Ab1=((a1.*Ac.*sind(si1-si))-(a1.*(Wc.^2).*cosd(si1-si))-(b1.*(Wb1.^2).*cosd(si1-beta1))-(c1.*(Wc1.^2)))./(b1.*sind(beta1-si1));
            Ac1=((a1.*Ac.*sind(beta1-si))-(a1.*(Wc.^2).*cosd(beta1-si))-(b1.*(Wb1.^2))+(c1.*(Wc1.^2).*cosd(beta1-si1)))./(c1.*sind(beta1-si1));
            
  
            %positions of each point
            O=[0;0];
            A=[a*cosd(theta);a*sind(theta)];
            B=[d*cosd(theta1)+c*cosd(si);d*sind(theta1)+c*sind(si)];
            C=[d*cosd(theta1);d*sind(theta1)];
            D=[d*cosd(theta1)+c*cosd(si)+b1*cosd(beta1);d*sind(theta1)+c*sind(si)+b1*sind(beta1)];
            E=[d*cosd(theta1)+d1*cosd(theta2);d*sind(theta1)+d1*sind(theta2)];
            
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
            B_v=[0 B_v]
            
            B_a=diff(B_v)./diff(time);
            B_a=[0 B_a];     
                      
                       
            D_x=D(1,:);
            D_y=D(2,:);
            
            D_vx=diff(D_x)./diff(time);
            D_vy=diff(D_y)./diff(time);
            
            D_v=sqrt(D_vx.^2+D_vy.^2);
            D_v=[0 D_v];
            
            D_a=diff(D_v)./diff(time);
            D_a=[0 D_a];
            
                
            as=app.asEditField.Value;   Fas=app.FasEditField.Value;     thetaFas=app.thetaFasEditField.Value;
            bs=app.bsEditField.Value;  Fbs=app.FbsEditField.Value;      thetaFbs=app.thetaFbsEditField.Value;
            cs=app.csEditField.Value;   Fcs=app.FcsEditField.Value;      thetaFcs=app.thetaFcsEditField.Value;
            a1s=cs;  Fa1s=Fcs;   thetaFa1s=thetaFcs;
            b1s=app.b1sEditField.Value;  Fb1s=app.Fb1sEditField.Value;     thetaFb1s=app.thetaFb1sEditField.Value;
            c1s=app.c1sEditField.Value; Fc1s=app.Fc1sEditField.Value;     thetaFc1s=app.thetaFc1sEditField.Value;  
        
            
            
            tor = Staticforce_A(app);
            
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
            d_1 = Document(fullFileName,"docx");
            open(d_1);
            if d_1 ~= -1
                %1st given data table
               tableStyle = { Width("110%"), ...
                               Border("solid"), ...
                               RowSep("solid"), ...
                               ColSep("solid") };

               append(d_1,Heading1("GIVEN DATA:"));
                BodyContent = {'Time range (secs)', time_of_running; ...
                               'Time step (secs)', delt; ...
                               'Initial angle made by crank1 with horizontal (degrees)',initial_angle; ...
                                'Angle made by fixed link1 with horizontal (degrees)' ,theta1; ... 
                                'Angle made by fixed link2 with horizontal (degrees)' ,theta2; ... 
                                'Anglular velocity of driver crank1 (rad/s) ' ,Wa};
                
                    tableContent_1 = [BodyContent];
                    
                    table = Table(tableContent_1);
                    table.Style = tableStyle;
                    
                    table.TableEntriesHAlign = "center";
                    append(d_1, table); 
                    
                    % 2nd - inputs table
                    tableStyle_1 = { Width("80%"), ...
                                   Border("solid"), ...
                                   RowSep("solid"), ...
                                   ColSep("solid") };
    
                    append(d_1,Heading1("INPUTS: "));
                    
                    HeaderContent = {'LINKS','LENGTH (metres)','TYPE'};
                    BodyContent = {'a',a,t1;'b',b,t2;'c',c,t3;...
                                    'd',d,t4;'b1',b1,t5;'c1',c1,t6;'d1',d1,t7};
                    tableContent_2 = [HeaderContent;BodyContent];
                    
                    table = Table(tableContent_2);
                    table.Style = tableStyle_1;
                    
                    table.TableEntriesHAlign = "center";
                    append(d_1, table);
                    
                     %3rd inputs table
                    append(d_1,Heading1("FORCES APPLIED ON LINKS"))
                    
                    formatSpec = 'On CRANK1 %d (N) at %d (m) at an angle %d (deg) with horizontal';                   
                    str2 = sprintf(formatSpec,Fas,as,thetaFas);
                    append(d_1,str2);
                    
                    formatSpec = 'On COUPLER1 %d (N) at %d (m) at an angle %d (deg) with horizontal';                   
                    str2 = sprintf(formatSpec,Fbs,bs,thetaFbs);
                    append(d_1,str2);
                    
                    formatSpec = 'On ROCKER1 %d (N) at %d (m) at an angle %d (deg) with horizontal';                   
                    str2 = sprintf(formatSpec,Fcs,cs,thetaFcs);
                    append(d_1,str2);
                    
                    formatSpec = 'On ROCKER2 %d (N) at %d (m) at an angle %d (deg) with horizontal';                   
                    str2 = sprintf(formatSpec,Fc1s,c1s,thetaFc1s);
                    append(d_1,str2);
                    
                    formatSpec = 'On COUPLER2 %d (N) at %d (m) at an angle %d (deg) with horizontal';                   
                    str2 = sprintf(formatSpec,Fb1s,b1s,thetaFb1s);
                    append(d_1,str2);                   
                   
                    
                    
                    %3rd - outputs table 1
                    
                    headerContent = {'Time (secs)','CRANK1 angle with horizontal (deg)',...
                        'ROCKER1 angle(deg)','COUPLER1 angle(deg)',...
                       'ROCKER2 angle(deg)','COUPLER2 angle(deg)'};
                    bodyContent = [time',theta',si',beta',si1',beta1'];
                    
                    data_str = string(bodyContent)
                    %round to 2 decimal places
                    for i = 1:numel(data_str)
                    data_str(i) = sprintf('%.2f',data_str(i))
                    end            
                    
                    tableContent = [headerContent; data_str];            
                    
                    append(d_1,Heading1("OUTPUTS - TABLE 1"));
                    
                    table = Table(tableContent);
                    table.Style = tableStyle;
                    
                    table.TableEntriesHAlign = "center";
                    append(d_1, table);
                                
            
                   % 3rd outputs table 2
                    table_Style = { Width("100%"), ...
                                   Border("solid"), ...
                                   RowSep("solid"), ...
                                   ColSep("solid") };

                    header_Content = {'Time (secs)','Velocity_A (m/s)','Velocity_B (m/s)','Velocity_D (m/s)',...
                              'Acceleration_A (m/s^2)','Acceleration_B (m/s^2)','Acceleration_D (m/s^2)' };
                    body_Content = [time',A_v',B_v',D_v',A_a',B_a',D_a'];                   
                    

                    
                    
                    data__str = string(body_Content);
                    %round to 2 decimal places
                    for i = 1:numel(data__str)
                    data__str(i) = sprintf('%.2f',data__str(i));
                    end
                    
                    
                    table_Content = [header_Content; data__str];
                    
                    
                    append(d_1,Heading1("OUTPUTS TABLE - 2"));
                    
                    table = Table(table_Content);
                    table.Style = table_Style;
                    
                    table.TableEntriesHAlign = "center";
                    append(d_1, table);
                    
                    formatSpec = 'Net torque acting on crank1 is %d N/m';
                    A1 = tor;
                    
                    str = sprintf(formatSpec,A1)
                    
                    append(d_1,Heading1(str));
                    
                    close(d_1);
                            
                
                
             else
                warningMessage = sprintf('Cannot open file:\n', fullFileName);
                uiwait(warndlg(warningMessage));
            end  
            
  
            
        end
        
        function [res] = Staticforce_A(app)
            
            l1 = app.LENGTHmEditField_8.Value;
            l2 = app.LENGTHmEditField_2.Value;
            l3= app.LENGTHmEditField_3.Value;
            l4= app.LENGTHmEditField_4.Value;
            
            l5 = app.LENGTHmEditField_5.Value;
            l6 = app.LENGTHmEditField_6.Value;
            l7 = app.LENGTHmEditField_7.Value;
            
            time_of_running = app.DurationofanimationsecsEditField.Value;
            delt =app.TimestepsecsEditField.Value;
            
            Wa = app.AngularvelocityofDrivercrankradsEditField.Value;
            theta1 = app.Inclinationoffixedlink1degreesEditField.Value;
            theta2 = app.Inclinationoffixedlink2degreesEditField.Value;
            initial_angle = app.Initialanglemadebycrank1withhorizontaldegreesEditField.Value;

            Aa=0;           
            
            final_angle=initial_angle+((time_of_running*Wa)*180/pi);
            
            theta=initial_angle:(Wa*180/pi)/2:final_angle;            
            time=linspace(0,time_of_running,length(theta));
            
            
            B=[l3;l5;l6;l7];
            Q=sort(B);
            
            A = [l1;l2;l3;l4];            
            P=sort(A);
            
            if(P(1)+P(4)<P(2)+P(3))
               a=P(1);
               b=P(4);
               c=P(2);
               d=P(3);
               if(Q(1)+Q(4)<Q(2)+Q(3))
                 if(c==Q(4))
                   a1=Q(4);
                   b1=Q(2);
                   c1=Q(3);
                   d1=Q(1); 
                 else 
                    b1=Q(1);        
                    c1=Q(4);
                    if(P(2)==Q(2))
                        a1=Q(2);
                        d1 = Q(3);
                    else d1 = Q(2);
                         a1 = Q(3);
                    end 
                 end
                 check_3=1;
               elseif(Q(1)+Q(4)>Q(2)+Q(3))
                   if(c==Q(1)) 
                        b1=Q(4)
                        c1=Q(2)
                        d1=Q(2)
                        a1=Q(1)
                        check_3=1;
                    else 
                      d1=Q(1);
                      a1=P(2);
                      b1=Q(4);
                      c1=Q(3);
                      check_3=1;
                   end  
                    
               else check_3 = 0;
                    
               end
            else check_3 =0;
                                 
            end
            
            as=app.asEditField.Value;   Fas=app.FasEditField.Value;     thetaFas=app.thetaFasEditField.Value;
            bs=app.bsEditField.Value;  Fbs=app.FbsEditField.Value;      thetaFbs=app.thetaFbsEditField.Value;
            cs=app.csEditField.Value;   Fcs=app.FcsEditField.Value;      thetaFcs=app.thetaFcsEditField.Value;
            a1s=cs;  Fa1s=Fcs;   thetaFa1s=thetaFcs;
            b1s=app.b1sEditField.Value;  Fb1s=app.Fb1sEditField.Value;     thetaFb1s=app.thetaFb1sEditField.Value;
            c1s=app.c1sEditField.Value; Fc1s=app.Fc1sEditField.Value;     thetaFc1s=app.thetaFc1sEditField.Value;  
            
            
            if as>a
                f1= errordlg("Analysis not possible");
                z=0;
            elseif bs>b
                f1 = errordlg("Analysis not possible");
                z=0;
            elseif cs>c
                f1= errordlg("Analysis not possible");
                z=0;
            elseif b1s>b1
                f1= errordlg("Analysis not possible");
                z=0;
             elseif c1s>c1
                f1= errordlg("Analysis not possible");
                z=0;
            else                
                z=1;
            end
            
            theta=app.Initialanglemadebycrank1withhorizontaldegreesEditField.Value; %defined by user
            
            if(z==1)
            k=((a^2-b^2+c^2+d^2)/2);
            A=-a.*(d.*cosd(theta1)-c).*cosd(theta)+k-(c*d*cosd(theta1))-a*d.*sind(theta).*sind(theta1);
            B=-(2.*a.*c.*sind(theta)-2*c*d.*sind(theta1));
            C=-a*(d*cosd(theta1)+c).*cosd(theta)+k+(c*d.*cosd(theta1))-a.*d.*sind(theta).*sind(theta1);
            si=2.*atand((-B-sqrt((B.^2)-4.*A.*C))./(2*A))
            
            G=((a^2+b^2-c^2+d^2)/2);
            D=-a.*(d.*cosd(theta1)+b).*cosd(theta)+G+(d*cosd(theta1)*b)-a*d*sind(theta1).*sind(theta);
            E=2*a*b.*sind(theta)-2*b*d.*sind(theta1);
            F=-(a.*(d*cosd(theta1)-b).*cosd(theta))+G-(b*d*cosd(theta1))-a.*d.*sind(theta).*sind(theta1);
            beta=2.*atand((-E+sqrt((E.^2)-4.*D.*F))./(2*D))
            
            k1=((a1^2-b1^2+c1^2+d1^2)/2);
            A1=-a1.*(d1.*cosd(theta2)-c1).*cosd(si)+k1-(c1*d1*cosd(theta2))-a1*d1.*sind(si).*sind(theta2);
            B1=-(2.*a1.*c1.*sind(si)-2*c1*d1.*sind(theta2));
            C1=-a1*(d1*cosd(theta2)+c1).*cosd(si)+k1+(c1*d1.*cosd(theta2))-a1.*d1.*sind(si).*sind(theta2);
            si1=2.*atand((-B1-sqrt((B1.^2)-4.*A1.*C1))./(2*A1))
            
            G1=((a1^2+b1^2-c1^2+d1^2)/2);
            D1=-a1.*(d1.*cosd(theta2)+b1).*cosd(si)+G1+(d1*cosd(theta2)*b1)-a1*d1*sind(theta2).*sind(si);
            E1=2*a1*b1.*sind(si)-2*b1*d1.*sind(theta2);
            F1=-(a1.*(d1*cosd(theta2)-b1).*cosd(si))+G1-(b1*d1*cosd(theta2))-a1.*d1.*sind(si).*sind(theta2);
            beta1=2.*atand((-E1+sqrt((E1.^2)-4.*D1.*F1))./(2*D1))
            
            
            syms Fab 
            moment_1=cross(as*[cosd(theta) sind(theta) 0],Fas*[cosd(thetaFas) sind(thetaFas) 0])+cross(a*[cosd(theta) sind(theta) 0],Fab*[cosd(beta) sind(beta) 0])==0
            moment_1 = moment_1(3);
            Fab = double(rhs((isolate(moment_1,Fab))))
            Torque_1 = cross( a*[cosd(theta) sind(theta) 0], Fab*[cosd(beta) sind(beta) 0])
            fprintf("%f Nm ",Torque_1(3))            
            
            
            syms Fab
            moment_2=cross(bs*[cosd(beta) sind(beta) 0],Fbs*[cosd(thetaFbs) sind(thetaFbs) 0])+cross(b*[cosd(beta) sind(beta) 0],Fab*[cosd(si) sind(si) 0])==0
            moment_2 = moment_2(3);
            Fab = double(rhs((isolate(moment_2,Fab))))
            Torque_2 = cross( a*[cosd(theta) sind(theta) 0], Fab*([cosd(si) sind(si) 0]))
            fprintf("%f Nm ",Torque_2(3))
            
            
            syms Fab 
            moment_3=cross(cs*[cosd(si) sind(si) 0],Fcs*[cosd(thetaFcs) sind(thetaFcs) 0])+cross([c*cosd(si) c*sind(si) 0],Fab*[cosd(beta) sind(beta) 0])==0
            moment_3 = moment_3(3);
            Fab = double(rhs((isolate(moment_3,Fab))))
            Torque_3 = cross( a*[cosd(theta) sind(theta) 0], Fab*([cosd(beta) sind(beta) 0]))
            fprintf("%f Nm ",Torque_3(3))
            
            syms Fab Fcb1
            moment_eq1=cross(bs*[cosd(beta1) sind(beta1) 0],Fb1s*[cosd(thetaFb1s) sind(thetaFb1s) 0])+cross(b1*[cosd(beta1) sind(beta1) 0],Fcb1*[cosd(si1) sind(si1) 0])==0
            moment_eq1 = moment_eq1(3);
            Fcb1 = double(rhs((isolate(moment_eq1,Fcb1))))
            moment_4=cross(Fcb1*[cosd(si1) sind(si1) 0],c*[cosd(si) sind(si) 0])+cross(Fab*[cosd(beta) sind(beta) 0],c*[cosd(si) sind(si) 0])==0
            moment_4 = moment_4(3);
            Fab = double(rhs((isolate(moment_4,Fab))))
            Torque_4 = cross( a*[cosd(theta) sind(theta) 0], Fab*([cosd(beta) sind(beta) 0]))
            fprintf("%f Nm ",Torque_4(3))
            
            syms Fab Fcb1
            moment_eq2=cross(cs*[cosd(si1) sind(si1) 0],Fc1s*[cosd(thetaFc1s) sind(thetaFc1s) 0])+cross(c1*[cosd(si) sind(si) 0],Fcb1*[cosd(beta1) sind(beta1) 0])==0
            moment_eq2 = moment_eq2(3);
            Fcb1 = double(rhs((isolate(moment_eq2,Fcb1))))
            moment_5=cross(Fcb1*[cosd(beta1) sind(beta1) 0],c*[cosd(si) sind(si) 0])+cross(Fab*[cosd(beta) sind(beta) 0],c*[cosd(si) sind(si) 0])==0
            moment_5 = moment_5(3);
            Fab = double(rhs((isolate(moment_5,Fab))))
            Torque_5 = cross( a*[cosd(theta) sind(theta) 0], Fab*([cosd(beta) sind(beta) 0]))
            fprintf("%f N-m ",Torque_5(3))
            
            NET_TORQUE(3)=Torque_1(3)+Torque_2(3)+Torque_3(3)+Torque_4(3)+Torque_5(3)
            fprintf("%f N-m ",NET_TORQUE(3))
            res=NET_TORQUE(3);
          
%           message = sprintf('Net torque acting on crank at equilibrium %d N/m',res);
%           uiwait(msgbox(message,'modal'));
            end
            
            return;
            
        end
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: CALCULATEButton
        function CALCULATEButtonPushed(app, event)
                        
           val = check_conditions_t1(app);
            if val 
                if app.VelocityandaccelerationanalysisCheckBox.Value
                     Velocity_accleration_A(app);
                elseif app.StaticforceanalysisCheckBox.Value                    
                     st=Staticforce_A(app);
                     message = sprintf('Net torque acting on crank at equilibrium %d N/m',st);
                     uiwait(msgbox(message,'modal'));
                 elseif app.AnimationCheckBox.Value
                     Animation_t1(app);
                end
            else f = errordlg('The given inputs do not the match the conditions required for the mechanism','Invalid Inputs');
                
            end
               
        end

        % Button pushed function: SAVEASButton
        function SAVEASButtonPushed(app, event)
            
           create_table(app);                    
           
            
        end

        % Button pushed function: RESETButton
        function RESETButtonPushed(app, event)
            %set(handles.my_edit_box,'String','');
%             clc;
%             clear all;
            close all;
            app.DurationofanimationsecsEditField.Value = 0;
            app.Initialanglemadebycrank1withhorizontaldegreesEditField.Value =0;
            app.Inclinationoffixedlink2degreesEditField.Value=0;
            app.AngularvelocityofDrivercrankradsEditField.Value=0;
            app.Inclinationoffixedlink1degreesEditField.Value=0;
            app.TimestepsecsEditField.Value=0;
            app.LENGTHmEditField_4.Value=0;
            app.LENGTHmEditField_3.Value=0;
            app.LENGTHmEditField_2.Value=0;
            app.LENGTHmEditField_8.Value=0;
            app.LENGTHmEditField_5.Value=0;
            app.LENGTHmEditField_6.Value=0;
            app.LENGTHmEditField_7.Value=0;
%             app.TYPEDropDown_13.Value = app.TYPEDropDown_13.Items(1);
%             app.TYPEDropDown_12.Value = app.TYPEDropDown_12.Items(1);
%             app.TYPEDropDown_11.Value = app.TYPEDropDown_11.Items(1);
%             app.TYPEDropDown_10.Value = app.TYPEDropDown_10.Items(1);
%             app.TYPEDropDown_9.Value = app.TYPEDropDown_9.Items(1);
%             app.TYPEDropDown_8.Value=app.TYPEDropDown_8.Items(1);
            app.NetaccelerationDropDown.Value=app.NetaccelerationDropDown.Items(1);
            app.NetvelocityDropDown.Value=app.NetvelocityDropDown.Items(1);
            app.AngularaccelarationDropDown.Value = app.AngularaccelarationDropDown.Items(1);
            app.AngularvelocityDropDown.Value=app.AngularvelocityDropDown.Items(1);
%             app.TYPEDropDown.Value=app.TYPEDropDown.Items(1);
            app.AnimationCheckBox.Value = false;
            app.VelocityandaccelerationanalysisCheckBox.Value = false;
            app.StaticforceanalysisCheckBox.Value = false;
            
        
        end

        % Button pushed function: CLOSEButton
        function CLOSEButtonPushed(app, event)
            delete(app);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 923 565];
            app.UIFigure.Name = 'MATLAB App';

            % Create INPUTLabel
            app.INPUTLabel = uilabel(app.UIFigure);
            app.INPUTLabel.Position = [190 524 106 26];
            app.INPUTLabel.Text = 'INPUT';

            % Create LENGTHmEditField_8Label
            app.LENGTHmEditField_8Label = uilabel(app.UIFigure);
            app.LENGTHmEditField_8Label.HorizontalAlignment = 'right';
            app.LENGTHmEditField_8Label.Position = [37 460 72 22];
            app.LENGTHmEditField_8Label.Text = 'LENGTH(m)';

            % Create LENGTHmEditField_8
            app.LENGTHmEditField_8 = uieditfield(app.UIFigure, 'numeric');
            app.LENGTHmEditField_8.Position = [118 456 69 29];

            % Create aLabel
            app.aLabel = uilabel(app.UIFigure);
            app.aLabel.Position = [21 456 25 30];
            app.aLabel.Text = 'a';

            % Create bLabel
            app.bLabel = uilabel(app.UIFigure);
            app.bLabel.Position = [21 410 25 23];
            app.bLabel.Text = 'b';

            % Create LENGTHmEditField_2Label
            app.LENGTHmEditField_2Label = uilabel(app.UIFigure);
            app.LENGTHmEditField_2Label.HorizontalAlignment = 'right';
            app.LENGTHmEditField_2Label.Position = [43 408 71 22];
            app.LENGTHmEditField_2Label.Text = 'LENGTH(m)';

            % Create LENGTHmEditField_2
            app.LENGTHmEditField_2 = uieditfield(app.UIFigure, 'numeric');
            app.LENGTHmEditField_2.Position = [123 404 69 29];

            % Create cLabel
            app.cLabel = uilabel(app.UIFigure);
            app.cLabel.Position = [20 368 25 23];
            app.cLabel.Text = 'c';

            % Create LENGTHmEditField_3Label
            app.LENGTHmEditField_3Label = uilabel(app.UIFigure);
            app.LENGTHmEditField_3Label.HorizontalAlignment = 'right';
            app.LENGTHmEditField_3Label.Position = [39 365 71 22];
            app.LENGTHmEditField_3Label.Text = 'LENGTH(m)';

            % Create LENGTHmEditField_3
            app.LENGTHmEditField_3 = uieditfield(app.UIFigure, 'numeric');
            app.LENGTHmEditField_3.Position = [119 362 69 29];

            % Create dLabel
            app.dLabel = uilabel(app.UIFigure);
            app.dLabel.Position = [18 320 25 23];
            app.dLabel.Text = 'd';

            % Create LENGTHmEditField_4Label
            app.LENGTHmEditField_4Label = uilabel(app.UIFigure);
            app.LENGTHmEditField_4Label.HorizontalAlignment = 'right';
            app.LENGTHmEditField_4Label.Position = [39 323 71 22];
            app.LENGTHmEditField_4Label.Text = 'LENGTH(m)';

            % Create LENGTHmEditField_4
            app.LENGTHmEditField_4 = uieditfield(app.UIFigure, 'numeric');
            app.LENGTHmEditField_4.Position = [119 320 69 29];

            % Create ANALYSISLabel
            app.ANALYSISLabel = uilabel(app.UIFigure);
            app.ANALYSISLabel.Position = [718 513 60 37];
            app.ANALYSISLabel.Text = 'ANALYSIS';

            % Create TIMEINTERVALLabel
            app.TIMEINTERVALLabel = uilabel(app.UIFigure);
            app.TIMEINTERVALLabel.Position = [331 485 96 40];
            app.TIMEINTERVALLabel.Text = 'TIME INTERVAL';

            % Create DurationofanimationsecsEditFieldLabel
            app.DurationofanimationsecsEditFieldLabel = uilabel(app.UIFigure);
            app.DurationofanimationsecsEditFieldLabel.HorizontalAlignment = 'right';
            app.DurationofanimationsecsEditFieldLabel.Position = [232 455 152 22];
            app.DurationofanimationsecsEditFieldLabel.Text = 'Duration of animation(secs)';

            % Create DurationofanimationsecsEditField
            app.DurationofanimationsecsEditField = uieditfield(app.UIFigure, 'numeric');
            app.DurationofanimationsecsEditField.Position = [399 451 48 30];

            % Create TimestepsecsEditFieldLabel
            app.TimestepsecsEditFieldLabel = uilabel(app.UIFigure);
            app.TimestepsecsEditFieldLabel.HorizontalAlignment = 'right';
            app.TimestepsecsEditFieldLabel.Position = [457 454 91 22];
            app.TimestepsecsEditFieldLabel.Text = 'Time step(secs)';

            % Create TimestepsecsEditField
            app.TimestepsecsEditField = uieditfield(app.UIFigure, 'numeric');
            app.TimestepsecsEditField.Position = [563 450 48 30];

            % Create ANGLEVARIATIONLabel
            app.ANGLEVARIATIONLabel = uilabel(app.UIFigure);
            app.ANGLEVARIATIONLabel.Position = [331 386 112 40];
            app.ANGLEVARIATIONLabel.Text = 'ANGLE VARIATION';

            % Create Inclinationoffixedlink1degreesEditFieldLabel
            app.Inclinationoffixedlink1degreesEditFieldLabel = uilabel(app.UIFigure);
            app.Inclinationoffixedlink1degreesEditFieldLabel.HorizontalAlignment = 'right';
            app.Inclinationoffixedlink1degreesEditFieldLabel.Position = [268 330 181 22];
            app.Inclinationoffixedlink1degreesEditFieldLabel.Text = 'Inclination of fixed link1(degrees)';

            % Create Inclinationoffixedlink1degreesEditField
            app.Inclinationoffixedlink1degreesEditField = uieditfield(app.UIFigure, 'numeric');
            app.Inclinationoffixedlink1degreesEditField.Position = [464 326 48 30];

            % Create VelocityandaccelerationanalysisCheckBox
            app.VelocityandaccelerationanalysisCheckBox = uicheckbox(app.UIFigure);
            app.VelocityandaccelerationanalysisCheckBox.Text = 'Velocity and acceleration analysis';
            app.VelocityandaccelerationanalysisCheckBox.Position = [666 450 203 32];

            % Create StaticforceanalysisCheckBox
            app.StaticforceanalysisCheckBox = uicheckbox(app.UIFigure);
            app.StaticforceanalysisCheckBox.Text = 'Static force analysis';
            app.StaticforceanalysisCheckBox.Position = [666 421 188 32];

            % Create CALCULATEButton
            app.CALCULATEButton = uibutton(app.UIFigure, 'push');
            app.CALCULATEButton.ButtonPushedFcn = createCallbackFcn(app, @CALCULATEButtonPushed, true);
            app.CALCULATEButton.Position = [522 19 137 31];
            app.CALCULATEButton.Text = 'CALCULATE';

            % Create SAVEASButton
            app.SAVEASButton = uibutton(app.UIFigure, 'push');
            app.SAVEASButton.ButtonPushedFcn = createCallbackFcn(app, @SAVEASButtonPushed, true);
            app.SAVEASButton.Position = [296 17 143 35];
            app.SAVEASButton.Text = 'SAVE AS';

            % Create RESETButton
            app.RESETButton = uibutton(app.UIFigure, 'push');
            app.RESETButton.ButtonPushedFcn = createCallbackFcn(app, @RESETButtonPushed, true);
            app.RESETButton.Position = [63 19 149 35];
            app.RESETButton.Text = 'RESET';

            % Create CLOSEButton
            app.CLOSEButton = uibutton(app.UIFigure, 'push');
            app.CLOSEButton.ButtonPushedFcn = createCallbackFcn(app, @CLOSEButtonPushed, true);
            app.CLOSEButton.Position = [723 21 142 35];
            app.CLOSEButton.Text = 'CLOSE';

            % Create AngularvelocityofDrivercrankradsEditFieldLabel
            app.AngularvelocityofDrivercrankradsEditFieldLabel = uilabel(app.UIFigure);
            app.AngularvelocityofDrivercrankradsEditFieldLabel.HorizontalAlignment = 'right';
            app.AngularvelocityofDrivercrankradsEditFieldLabel.Position = [263 267 204 29];
            app.AngularvelocityofDrivercrankradsEditFieldLabel.Text = 'Angular velocity of Driver crank(rad/s)';

            % Create AngularvelocityofDrivercrankradsEditField
            app.AngularvelocityofDrivercrankradsEditField = uieditfield(app.UIFigure, 'numeric');
            app.AngularvelocityofDrivercrankradsEditField.Position = [466 268 45 30];

            % Create b1Label
            app.b1Label = uilabel(app.UIFigure);
            app.b1Label.Position = [18 285 25 23];
            app.b1Label.Text = 'b1';

            % Create c1Label
            app.c1Label = uilabel(app.UIFigure);
            app.c1Label.Position = [15 245 25 23];
            app.c1Label.Text = 'c1';

            % Create d1Label
            app.d1Label = uilabel(app.UIFigure);
            app.d1Label.Position = [15 206 58 23];
            app.d1Label.Text = 'd1';

            % Create LENGTHmEditField_5Label
            app.LENGTHmEditField_5Label = uilabel(app.UIFigure);
            app.LENGTHmEditField_5Label.HorizontalAlignment = 'right';
            app.LENGTHmEditField_5Label.Position = [43 285 71 22];
            app.LENGTHmEditField_5Label.Text = 'LENGTH(m)';

            % Create LENGTHmEditField_5
            app.LENGTHmEditField_5 = uieditfield(app.UIFigure, 'numeric');
            app.LENGTHmEditField_5.Position = [123 282 69 29];

            % Create LENGTHmEditField_6Label
            app.LENGTHmEditField_6Label = uilabel(app.UIFigure);
            app.LENGTHmEditField_6Label.HorizontalAlignment = 'right';
            app.LENGTHmEditField_6Label.Position = [38 242 71 22];
            app.LENGTHmEditField_6Label.Text = 'LENGTH(m)';

            % Create LENGTHmEditField_6
            app.LENGTHmEditField_6 = uieditfield(app.UIFigure, 'numeric');
            app.LENGTHmEditField_6.Position = [118 239 69 29];

            % Create LENGTHmEditField_7Label
            app.LENGTHmEditField_7Label = uilabel(app.UIFigure);
            app.LENGTHmEditField_7Label.HorizontalAlignment = 'right';
            app.LENGTHmEditField_7Label.Position = [39 203 71 22];
            app.LENGTHmEditField_7Label.Text = 'LENGTH(m)';

            % Create LENGTHmEditField_7
            app.LENGTHmEditField_7 = uieditfield(app.UIFigure, 'numeric');
            app.LENGTHmEditField_7.Position = [119 200 69 29];

            % Create Inclinationoffixedlink2degreesEditFieldLabel
            app.Inclinationoffixedlink2degreesEditFieldLabel = uilabel(app.UIFigure);
            app.Inclinationoffixedlink2degreesEditFieldLabel.HorizontalAlignment = 'right';
            app.Inclinationoffixedlink2degreesEditFieldLabel.Position = [268 301 181 22];
            app.Inclinationoffixedlink2degreesEditFieldLabel.Text = 'Inclination of fixed link2(degrees)';

            % Create Inclinationoffixedlink2degreesEditField
            app.Inclinationoffixedlink2degreesEditField = uieditfield(app.UIFigure, 'numeric');
            app.Inclinationoffixedlink2degreesEditField.Position = [464 297 48 30];

            % Create AnimationCheckBox
            app.AnimationCheckBox = uicheckbox(app.UIFigure);
            app.AnimationCheckBox.Text = 'Animation';
            app.AnimationCheckBox.Position = [666 390 188 32];

            % Create Initialanglemadebycrank1withhorizontaldegreesEditFieldLabel
            app.Initialanglemadebycrank1withhorizontaldegreesEditFieldLabel = uilabel(app.UIFigure);
            app.Initialanglemadebycrank1withhorizontaldegreesEditFieldLabel.HorizontalAlignment = 'right';
            app.Initialanglemadebycrank1withhorizontaldegreesEditFieldLabel.Position = [232 359 298 22];
            app.Initialanglemadebycrank1withhorizontaldegreesEditFieldLabel.Text = 'Initial angle made by crank1 with horizontal(degrees)';

            % Create Initialanglemadebycrank1withhorizontaldegreesEditField
            app.Initialanglemadebycrank1withhorizontaldegreesEditField = uieditfield(app.UIFigure, 'numeric');
            app.Initialanglemadebycrank1withhorizontaldegreesEditField.Position = [537 355 48 30];

            % Create PLOTSLabel
            app.PLOTSLabel = uilabel(app.UIFigure);
            app.PLOTSLabel.Position = [723 360 73 24];
            app.PLOTSLabel.Text = 'PLOTS';

            % Create AngularvelocityDropDownLabel
            app.AngularvelocityDropDownLabel = uilabel(app.UIFigure);
            app.AngularvelocityDropDownLabel.HorizontalAlignment = 'right';
            app.AngularvelocityDropDownLabel.Position = [658 334 90 22];
            app.AngularvelocityDropDownLabel.Text = 'Angular velocity';

            % Create AngularvelocityDropDown
            app.AngularvelocityDropDown = uidropdown(app.UIFigure);
            app.AngularvelocityDropDown.Items = {'Select', 'ROCKER1', 'COUPLER1', 'ROCKER2', 'COUPLER2'};
            app.AngularvelocityDropDown.Position = [763 334 106 22];
            app.AngularvelocityDropDown.Value = 'Select';

            % Create AngularaccelarationDropDownLabel
            app.AngularaccelarationDropDownLabel = uilabel(app.UIFigure);
            app.AngularaccelarationDropDownLabel.HorizontalAlignment = 'right';
            app.AngularaccelarationDropDownLabel.Position = [633 313 115 22];
            app.AngularaccelarationDropDownLabel.Text = 'Angular accelaration';

            % Create AngularaccelarationDropDown
            app.AngularaccelarationDropDown = uidropdown(app.UIFigure);
            app.AngularaccelarationDropDown.Items = {'Select', 'ROCKER1', 'COUPLER1', 'ROCKER2', 'COUPLER2'};
            app.AngularaccelarationDropDown.Position = [763 313 106 22];
            app.AngularaccelarationDropDown.Value = 'Select';

            % Create NetvelocityDropDownLabel
            app.NetvelocityDropDownLabel = uilabel(app.UIFigure);
            app.NetvelocityDropDownLabel.HorizontalAlignment = 'right';
            app.NetvelocityDropDownLabel.Position = [681 289 67 22];
            app.NetvelocityDropDownLabel.Text = 'Net velocity';

            % Create NetvelocityDropDown
            app.NetvelocityDropDown = uidropdown(app.UIFigure);
            app.NetvelocityDropDown.Items = {'Select', 'Point A', 'Point B', 'Point D'};
            app.NetvelocityDropDown.Position = [763 289 106 22];
            app.NetvelocityDropDown.Value = 'Select';

            % Create NetaccelerationDropDownLabel
            app.NetaccelerationDropDownLabel = uilabel(app.UIFigure);
            app.NetaccelerationDropDownLabel.HorizontalAlignment = 'right';
            app.NetaccelerationDropDownLabel.Position = [656 268 92 22];
            app.NetaccelerationDropDownLabel.Text = 'Net acceleration';

            % Create NetaccelerationDropDown
            app.NetaccelerationDropDown = uidropdown(app.UIFigure);
            app.NetaccelerationDropDown.Items = {'Select', 'Point A', 'PointB', 'Point D'};
            app.NetaccelerationDropDown.Position = [763 268 106 22];
            app.NetaccelerationDropDown.Value = 'Select';

            % Create FasEditFieldLabel
            app.FasEditFieldLabel = uilabel(app.UIFigure);
            app.FasEditFieldLabel.HorizontalAlignment = 'right';
            app.FasEditFieldLabel.Position = [311 221 26 22];
            app.FasEditFieldLabel.Text = 'Fas';

            % Create FasEditField
            app.FasEditField = uieditfield(app.UIFigure, 'numeric');
            app.FasEditField.Position = [349 216 35 32];

            % Create FbsEditFieldLabel
            app.FbsEditFieldLabel = uilabel(app.UIFigure);
            app.FbsEditFieldLabel.HorizontalAlignment = 'right';
            app.FbsEditFieldLabel.Position = [311 190 26 22];
            app.FbsEditFieldLabel.Text = 'Fbs';

            % Create FbsEditField
            app.FbsEditField = uieditfield(app.UIFigure, 'numeric');
            app.FbsEditField.Position = [349 185 35 32];

            % Create FcsEditFieldLabel
            app.FcsEditFieldLabel = uilabel(app.UIFigure);
            app.FcsEditFieldLabel.HorizontalAlignment = 'right';
            app.FcsEditFieldLabel.Position = [312 159 25 22];
            app.FcsEditFieldLabel.Text = 'Fcs';

            % Create FcsEditField
            app.FcsEditField = uieditfield(app.UIFigure, 'numeric');
            app.FcsEditField.Position = [349 154 35 32];

            % Create Fb1sEditFieldLabel
            app.Fb1sEditFieldLabel = uilabel(app.UIFigure);
            app.Fb1sEditFieldLabel.HorizontalAlignment = 'right';
            app.Fb1sEditFieldLabel.Position = [306 127 32 22];
            app.Fb1sEditFieldLabel.Text = 'Fb1s';

            % Create Fb1sEditField
            app.Fb1sEditField = uieditfield(app.UIFigure, 'numeric');
            app.Fb1sEditField.Position = [350 122 35 32];

            % Create Fc1sEditFieldLabel
            app.Fc1sEditFieldLabel = uilabel(app.UIFigure);
            app.Fc1sEditFieldLabel.HorizontalAlignment = 'right';
            app.Fc1sEditFieldLabel.Position = [307 96 31 22];
            app.Fc1sEditFieldLabel.Text = 'Fc1s';

            % Create Fc1sEditField
            app.Fc1sEditField = uieditfield(app.UIFigure, 'numeric');
            app.Fc1sEditField.Position = [350 91 35 32];

            % Create asEditFieldLabel
            app.asEditFieldLabel = uilabel(app.UIFigure);
            app.asEditFieldLabel.HorizontalAlignment = 'right';
            app.asEditFieldLabel.Position = [406 221 25 22];
            app.asEditFieldLabel.Text = 'as';

            % Create asEditField
            app.asEditField = uieditfield(app.UIFigure, 'numeric');
            app.asEditField.Position = [443 216 35 32];

            % Create bsEditFieldLabel
            app.bsEditFieldLabel = uilabel(app.UIFigure);
            app.bsEditFieldLabel.HorizontalAlignment = 'right';
            app.bsEditFieldLabel.Position = [405 190 25 22];
            app.bsEditFieldLabel.Text = 'bs';

            % Create bsEditField
            app.bsEditField = uieditfield(app.UIFigure, 'numeric');
            app.bsEditField.Position = [442 185 35 32];

            % Create csEditFieldLabel
            app.csEditFieldLabel = uilabel(app.UIFigure);
            app.csEditFieldLabel.HorizontalAlignment = 'right';
            app.csEditFieldLabel.Position = [405 159 25 22];
            app.csEditFieldLabel.Text = 'cs';

            % Create csEditField
            app.csEditField = uieditfield(app.UIFigure, 'numeric');
            app.csEditField.Position = [442 154 35 32];

            % Create b1sEditFieldLabel
            app.b1sEditFieldLabel = uilabel(app.UIFigure);
            app.b1sEditFieldLabel.HorizontalAlignment = 'right';
            app.b1sEditFieldLabel.Position = [405 127 25 22];
            app.b1sEditFieldLabel.Text = 'b1s';

            % Create b1sEditField
            app.b1sEditField = uieditfield(app.UIFigure, 'numeric');
            app.b1sEditField.Position = [442 122 35 32];

            % Create c1sEditFieldLabel
            app.c1sEditFieldLabel = uilabel(app.UIFigure);
            app.c1sEditFieldLabel.HorizontalAlignment = 'right';
            app.c1sEditFieldLabel.Position = [405 96 25 22];
            app.c1sEditFieldLabel.Text = 'c1s';

            % Create c1sEditField
            app.c1sEditField = uieditfield(app.UIFigure, 'numeric');
            app.c1sEditField.Position = [442 91 35 32];

            % Create thetaFasEditFieldLabel
            app.thetaFasEditFieldLabel = uilabel(app.UIFigure);
            app.thetaFasEditFieldLabel.HorizontalAlignment = 'right';
            app.thetaFasEditFieldLabel.Position = [494 221 52 22];
            app.thetaFasEditFieldLabel.Text = 'thetaFas';

            % Create thetaFasEditField
            app.thetaFasEditField = uieditfield(app.UIFigure, 'numeric');
            app.thetaFasEditField.Position = [558 216 35 32];

            % Create thetaFbsEditFieldLabel
            app.thetaFbsEditFieldLabel = uilabel(app.UIFigure);
            app.thetaFbsEditFieldLabel.HorizontalAlignment = 'right';
            app.thetaFbsEditFieldLabel.Position = [494 190 52 22];
            app.thetaFbsEditFieldLabel.Text = 'thetaFbs';

            % Create thetaFbsEditField
            app.thetaFbsEditField = uieditfield(app.UIFigure, 'numeric');
            app.thetaFbsEditField.Position = [558 185 35 32];

            % Create thetaFcsEditFieldLabel
            app.thetaFcsEditFieldLabel = uilabel(app.UIFigure);
            app.thetaFcsEditFieldLabel.HorizontalAlignment = 'right';
            app.thetaFcsEditFieldLabel.Position = [495 159 51 22];
            app.thetaFcsEditFieldLabel.Text = 'thetaFcs';

            % Create thetaFcsEditField
            app.thetaFcsEditField = uieditfield(app.UIFigure, 'numeric');
            app.thetaFcsEditField.Position = [558 154 35 32];

            % Create thetaFb1sEditFieldLabel
            app.thetaFb1sEditFieldLabel = uilabel(app.UIFigure);
            app.thetaFb1sEditFieldLabel.HorizontalAlignment = 'right';
            app.thetaFb1sEditFieldLabel.Position = [487 128 59 22];
            app.thetaFb1sEditFieldLabel.Text = 'thetaFb1s';

            % Create thetaFb1sEditField
            app.thetaFb1sEditField = uieditfield(app.UIFigure, 'numeric');
            app.thetaFb1sEditField.Position = [558 123 35 32];

            % Create thetaFc1sEditFieldLabel
            app.thetaFc1sEditFieldLabel = uilabel(app.UIFigure);
            app.thetaFc1sEditFieldLabel.HorizontalAlignment = 'right';
            app.thetaFc1sEditFieldLabel.Position = [487 97 58 22];
            app.thetaFc1sEditFieldLabel.Text = 'thetaFc1s';

            % Create thetaFc1sEditField
            app.thetaFc1sEditField = uieditfield(app.UIFigure, 'numeric');
            app.thetaFc1sEditField.Position = [557 92 35 32];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Variable_Stroke_mechanism_t1_exported

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
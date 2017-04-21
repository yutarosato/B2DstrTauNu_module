#ifndef USERINFO_H
#define USERINFO_H

#include "belle.h"
#include "particle/ParticleUserInfo.h"

#include "CLHEP/Vector/LorentzVector.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
  
  class UserInfo : public ParticleUserInfo
  {
  public:
    UserInfo();                                // Default constructor
    UserInfo( const  UserInfo& );              // Copy constructor 
    virtual ~UserInfo();                       // Destructor
    UserInfo* clone(void) const;               // constructs self object.
    UserInfo& operator = ( const UserInfo& );  // Copy operator


    // mode
    void rec_mode      ( const int&              v ){ m_rec_mode=v;       }
    void cntid         ( const int&              v ){ m_cntid=v;          }
    // particle
    void self          ( const int&              v ){ m_self=v;           }
    void m_org         ( const double&           v ){ m_m_org=v;          }
    void dm            ( const double&           v ){ m_dm=v;             }
    void m_pse         ( const double&           v ){ m_m_pse=v;          }
    void id            ( const int&              v ){ m_id=v;             }
    void motherid      ( const int&              v ){ m_motherid=v;       }
    void flavor        ( const int&              v ){ m_flavor=v;         }
    void angle         ( const double&           v ){ m_angle=v;          }
    void cos           ( const double&           v ){ m_cos=v;            }
    void Vcm           ( const HepLorentzVector& v ){ m_Vcm=v;            }
    void dr            ( const double&           v ){ m_dr=v;             }
    void dz            ( const double&           v ){ m_dz=v;             }
    void eidProb       ( const double&           v ){ m_eidProb=v;        }
    void muonLikelihood( const double&           v ){ m_muonLikelihood=v; }
    void selKPI        ( const double&           v ){ m_selKPI=v;         }
    void radgam        ( const int&              v ){ m_radgam=v;         }
    void radgam_n      ( const int&              v ){ m_radgam_n=v;       }
    // kinematic fit
    void kf            ( const int&              v ){ m_kf=v;             }
    void kf_cl         ( const double&           v ){ m_kf_cl=v;          }
    void kf_chisq      ( const double&           v ){ m_kf_chisq=v;       }
    void kf_ndf        ( const int&              v ){ m_kf_ndf=v;         }
   // B
    void thrust_angle  ( const double&           v ){ m_thrust_angle=v;   }
    void R2            ( const double&           v ){ m_R2=v;             }
    void Mbc           ( const double&           v ){ m_Mbc=v;            }
    void dE            ( const double&           v ){ m_dE=v;             }
    // tmp
    void tmp_int       ( const int&              v ){ m_tmp_int=v;        }
    void tmp_double    ( const double&           v ){ m_tmp_double=v;     }


    const int&              rec_mode(void)       const { return m_rec_mode;       }
    const int&              cntid(void)          const { return m_cntid;          }
    const int&              self(void)           const { return m_self;           }
    const double&           m_org(void)          const { return m_m_org;          }
    const double&           dm(void)             const { return m_dm;             }
    const double&           m_pse(void)          const { return m_m_pse;          }
    const int&              id(void)             const { return m_id;             }
    const int&              motherid(void)       const { return m_motherid;       }
    const int&              flavor(void)         const { return m_flavor;         }
    const double&           angle(void)          const { return m_angle;          }
    const double&           cos(void)            const { return m_cos;            }
    const HepLorentzVector& Vcm(void)            const { return m_Vcm;            }
    const double&           dr(void)             const { return m_dr;             }
    const double&           dz(void)             const { return m_dz;             }
    const double&           eidProb(void)        const { return m_eidProb;        }
    const double&           muonLikelihood(void) const { return m_muonLikelihood; }
    const double&           selKPI(void)         const { return m_selKPI;         }
    const int&              radgam(void)         const { return m_radgam;         }
    const int&              radgam_n(void)       const { return m_radgam_n;       }
    const int&              kf(void)             const { return m_kf;             }
    const double&           kf_cl(void)          const { return m_kf_cl;          }
    const double&           kf_chisq(void)       const { return m_kf_chisq;       }
    const int&              kf_ndf(void)         const { return m_kf_ndf;         }
    const double&           thrust_angle(void)   const { return m_thrust_angle;   }
    const double&           R2(void)             const { return m_R2;             }
    const double&           Mbc(void)            const { return m_Mbc;            }
    const double&           dE(void)             const { return m_dE;             }
    const int&              tmp_int(void)        const { return m_tmp_int;        }
    const double&           tmp_double(void)     const { return m_tmp_double;     }
    
  private:
    int              m_rec_mode;
    int              m_cntid;
    int              m_self;
    double           m_m_org;
    double           m_dm;
    double           m_m_pse;
    int              m_id;
    int              m_motherid;
    int              m_flavor;
    double           m_angle;
    double           m_cos;
    HepLorentzVector m_Vcm;
    double           m_dr;
    double           m_dz;
    double           m_eidProb;
    double           m_muonLikelihood;
    double           m_selKPI;
    int              m_radgam;
    int              m_radgam_n;
    int              m_kf;
    double           m_kf_cl;
    double           m_kf_chisq;
    int              m_kf_ndf;
    double           m_thrust_angle;
    double           m_R2;
    double           m_Mbc;
    double           m_dE;
    int              m_tmp_int;
    double           m_tmp_double;
  };
  
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
  


#include "UserInfo.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  // Default constructor
  UserInfo::UserInfo() :
    m_rec_mode(-1),
    m_cntid(-1),
    m_self(-1),
    m_m_org(-1),
    m_dm(-1),
    m_m_pse(-1),
    m_id(0),
    m_motherid(0),
    m_flavor(0),
    m_angle(-10),
    m_cos(-10),
    m_cos_ctrl(-10),
    m_Vcm(0,0,0,0),
    m_dr(-1.0),
    m_dz(-1.0),
    m_eidProb(2.0),
    m_muonLikelihood(2.0),
    m_selKPI(2.0),
    m_radgam(-1),
    m_radgam_n(0),
    m_kf(-1),
    m_kf_cl(-1.0),
    m_kf_chisq(-1.0),
    m_kf_ndf(-1),
    m_thrust_angle(-1.0),
    m_R2(-1.0),
    m_Mbc(-10),
    m_dE(-10),
    m_tmp_int(-10),
    m_tmp_double(-10)
  {}
  
  // Destructor
  UserInfo::~UserInfo(){}

  // Copy constructor 
  UserInfo::UserInfo(const UserInfo& x) :
    m_rec_mode(x.m_rec_mode),
    m_cntid(x.m_cntid),
    m_self(x.m_self),
    m_m_org(x.m_m_org),
    m_dm(x.m_dm),
    m_m_pse(x.m_m_pse),
    m_id(x.m_id),
    m_motherid(x.m_motherid),
    m_flavor(x.m_flavor),
    m_angle(x.m_angle),
    m_cos(x.m_cos),
    m_cos_ctrl(x.m_cos_ctrl),
    m_Vcm(x.m_Vcm),
    m_dr(x.m_dr),
    m_dz(x.m_dz),
    m_eidProb(x.m_eidProb),
    m_muonLikelihood(x.m_muonLikelihood),
    m_selKPI(x.m_selKPI),
    m_radgam(x.m_radgam),
    m_radgam_n(x.m_radgam_n),
    m_kf(x.m_kf),
    m_kf_cl(x.m_kf_cl),
    m_kf_chisq(x.m_kf_chisq),
    m_kf_ndf(x.m_kf_ndf),
    m_thrust_angle(x.m_thrust_angle),
    m_R2(x.m_R2),
    m_Mbc(x.m_Mbc),
    m_dE(x.m_dE),
    m_tmp_int(x.m_tmp_int),
    m_tmp_double(x.m_tmp_double)
  {}
  
  // constructs self object.
  UserInfo* UserInfo::clone(void) const
  {
    return (new UserInfo(*this));
  }

  // Copy operator
  UserInfo& UserInfo::operator = (const UserInfo& x)
  {
    ParticleUserInfo::operator = (x);
    m_rec_mode       = x.m_rec_mode;
    m_cntid          = x.m_cntid;
    m_self           = x.m_self;
    m_m_org          = x.m_m_org;
    m_dm             = x.m_dm;
    m_m_pse          = x.m_m_pse;
    m_id             = x.m_id;
    m_motherid       = x.m_motherid;
    m_flavor         = x.m_flavor;
    m_angle          = x.m_angle;
    m_cos            = x.m_cos;
    m_cos_ctrl       = x.m_cos_ctrl;
    m_Vcm            = x.m_Vcm;
    m_dr             = x.m_dr;
    m_dz             = x.m_dz;
    m_eidProb        = x.m_eidProb;
    m_muonLikelihood = x.m_muonLikelihood;
    m_selKPI         = x.m_selKPI;
    m_radgam         = x.m_radgam;
    m_radgam_n       = x.m_radgam_n;
    m_kf             = x.m_kf;
    m_kf_cl          = x.m_kf_cl;
    m_kf_chisq       = x.m_kf_chisq;
    m_kf_ndf         = x.m_kf_ndf;
    m_thrust_angle   = x.m_thrust_angle;
    m_R2             = x.m_R2;
    m_Mbc            = x.m_Mbc;
    m_dE             = x.m_dE;
    m_tmp_int        = x.m_tmp_int;
    m_tmp_double     = x.m_tmp_double;
  }
  
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

  

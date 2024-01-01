import React from "react";
import { BsChevronDoubleLeft, BsChevronDoubleRight, BsGithub, BsBugFill } from "react-icons/bs";

class OpenPageButton extends React.Component {
    openNewTab = () => { window.open(this.props.url, "_blank"); };
  
    render() {
        return (
            <button className={`square-button ${this.props.className}`} onClick={this.openNewTab}>
                {this.props.icon}
                <span>{this.props.title}</span>
            </button>
        );
    };
};

const Sidebar = props => {
    const sidebarOpen = <BsChevronDoubleLeft className="sidebar-toggle-icon" />;
    const sidebarClosed = <BsChevronDoubleRight className="sidebar-toggle-icon" />;

    return (
        <div className={"sidebar" + (props.isOpen ? " open" : "")}>
            <button onClick={props.toggleSidebar} className="sidebar-toggle">
                {props.isOpen ? sidebarOpen : sidebarClosed }
            </button>
            
            <div className="content">
                <div className="version">Version 0.1.0</div>
                <OpenPageButton 
                    className="sidebar-button" 
                    icon={<BsGithub />} 
                    title={"Open GitHub page"} 
                    url="https://github.com/moltools/CineMol" 
                />
                <OpenPageButton 
                    className="sidebar-button" 
                    icon={<BsBugFill />} 
                    title={"File bug report or request feature"} 
                    url="https://github.com/moltools/CineMol/issues" 
                />
            </div>
        </div>
    );
};

export default Sidebar;
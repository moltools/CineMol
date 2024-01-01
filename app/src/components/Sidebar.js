import React from "react";

const Sidebar = props => {
    const sidebarClass = props.isOpen ? "sidebar open" : "sidebar";

    return (
        <div className={sidebarClass}>
            <button onClick={props.toggleSidebar} className="sidebar-toggle">
                Toggle Sidebar
            </button>
        </div>
    );
};

export default Sidebar;